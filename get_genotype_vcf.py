#!/usr/bin/env python3

"""
This script will

- extract common variants
- extract rare variants
- export (both) as VCF files (one per chrom)
- also create a plink file for a subset of variants
for variance ratio estimation

this will be used as input for the
SAIGE-QTL association pipeline.

analysis-runner \
   --description "get common and rare variant VCFs" \
   --dataset "tenk10k" \
   --access-level "standard" \
   --output-dir saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/genotypes/ \
    python3 get_genotype_vcf.py --vds-path=gs://cpg-tenk10k-main/vds/tenk10k-genome-2-0.vds --chromosomes chr2 \
    --relateds-to-drop-path=gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/relateds_to_drop.ht \
    --qc-pass-variants-path=gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/variants_qc.ht
"""

import logging

import click

import hail as hl

from hail.methods import export_vcf

from cpg_utils import to_path
from cpg_utils.config import get_config, output_path
from cpg_utils.hail_batch import get_batch, init_batch


def can_reuse(path: str):
    """
    checks for existence of a Path
    if the path is a MT or VDS, checks for the success file
    """

    if not path:
        return False

    path_as_path = to_path(path)

    if path_as_path.suffix in ['.mt', '.ht']:
        path_as_path /= '_SUCCESS'
    if path_as_path.suffix in ['.vds']:
        path_as_path /= 'variant_data/_SUCCESS'

    return path_as_path.exists()


def add_remove_chr_and_index_job(vcf_path):
    """
    Reads a VCF file, it creates an .csi index file and
    and removes "chr" from the chromosome values

    Args:
    vcf_path: input path
    """
    # remove chr & add index file using bcftools
    vcf_input = get_batch().read_input(vcf_path)

    vcf_size = to_path(vcf_path).stat().st_size
    storage_required = ((vcf_size // 1024**3) or 1) * 2.2
    bcftools_job = get_batch().new_job(name='remove chr and index vcf')
    bcftools_job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.csi': '{root}.vcf.bgz.csi',
        }
    )
    bcftools_job.image(get_config()['images']['bcftools'])
    bcftools_job.cpu(4)
    bcftools_job.storage(f'{storage_required}Gi')
    # now remove "chr" from chromosome names using bcftools
    bcftools_job.command(
        'for num in {1..22} X Y; do echo "chr${num} ${num}" >> chr_update.txt; done'
    )
    bcftools_job.command(
        f"""
        bcftools annotate --rename-chrs chr_update.txt --set-id +'%CHROM\:%POS\:%REF\:%FIRST_ALT' {vcf_input} | \\
        bgzip -c > {bcftools_job.output['vcf.bgz']}
        bcftools index -c {bcftools_job.output['vcf.bgz']}
    """
    )
    logging.info('VCF rename/index jobs scheduled!')

    # save both output files
    get_batch().write_output(bcftools_job.output, vcf_path.removesuffix('.vcf.bgz'))


# inputs:
@click.option('--vds-path', help=' GCP gs:// path to the VDS')
@click.option('--chromosomes', help=' e.g., chr22,chrX ')
@click.option(
    '--relateds-to-drop-path',
    default='gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/relateds_to_drop.ht',
)
@click.option(
    '--qc-pass-variants-path',
    default='gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/variants_qc.ht',
)
@click.option('--cv-maf-threshold', default=0.01)
@click.option('--rv-maf-threshold', default=0.01)
@click.option('--exclude-multiallelic', is_flag=False)
@click.option('--exclude-indels', is_flag=False)
@click.command()
def main(
    vds_path: str,
    chromosomes: str,
    relateds_to_drop_path: str,
    qc_pass_variants_path: str,
    cv_maf_threshold: float,
    rv_maf_threshold: float,
    exclude_multiallelic: bool,
    exclude_indels: bool,
):
    """
    Write chromosome-level genotypes as VCF
    """

    init_batch(worker_memory='highmem', driver_memory='highmem')

    vds = hl.vds.read_vds(vds_path)
    vds_name = vds_path.split('/')[-1].split('.')[0]

    for chromosome in chromosomes.split(','):
        # create paths and check if they exist already
        cv_vcf_path = output_path(
            f'vds-{vds_name}/{chromosome}_common_variants.vcf.bgz'
        )
        cv_vcf_existence_outcome = can_reuse(cv_vcf_path)
        logging.info(f'Does {cv_vcf_path} exist? {cv_vcf_existence_outcome}')

        cv_mt_path = output_path(f'vds-{vds_name}/{chromosome}_common_variants.mt')
        cv_mt_existence_outcome = can_reuse(cv_mt_path)
        logging.info(f'Does {cv_mt_path} exist? {cv_mt_existence_outcome}')

        rv_vcf_path = output_path(f'vds-{vds_name}/{chromosome}_rare_variants.vcf.bgz')
        rv_vcf_existence_outcome = can_reuse(rv_vcf_path)
        logging.info(f'Does {rv_vcf_path} exist? {rv_vcf_existence_outcome}')

        rv_mt_path = output_path(f'vds-{vds_name}/{chromosome}_rare_variants.mt')
        rv_mt_existence_outcome = can_reuse(rv_mt_path)
        logging.info(f'Does {rv_mt_path} exist? {rv_mt_existence_outcome}')

        if (
            not cv_vcf_existence_outcome
            or not cv_mt_existence_outcome
            or not rv_vcf_existence_outcome
            or not rv_mt_existence_outcome
        ):
            # consider only relevant chromosome
            chrom_vds = hl.vds.filter_chromosomes(vds, keep=chromosome)

            # split multiallelic loci (necessary pre-densifying)
            chrom_vds = hl.vds.split_multi(chrom_vds, filter_changed_loci=True)

            # densify to matrix table object
            mt = hl.vds.to_dense_mt(chrom_vds)

            # filter out related samples
            related_ht = hl.read_table(relateds_to_drop_path)
            related_samples = hl.literal(related_ht.s.collect())
            mt = mt.filter_cols(~related_samples.contains(mt['s']))

            # filter out loci & variant QC
            qc_pass_variants = hl.read_table(qc_pass_variants_path)
            mt = mt.semi_join_rows(qc_pass_variants)
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # remove hom-ref
            if exclude_multiallelic:  # biallelic only (exclude multiallelic)
                mt = mt.filter_rows(~(mt.was_split))
            if exclude_indels:  # SNPs only (exclude indels)
                mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

            mt = hl.variant_qc(mt)

            if not rv_vcf_existence_outcome or not rv_mt_existence_outcome:
                # rare variants only
                rv_mt = mt.filter_rows(hl.min(mt.variant_qc.AF) < rv_maf_threshold)

                if not rv_mt_existence_outcome:
                    # save chrom + rare variant mt for group file script
                    rv_mt.write(rv_mt_path, overwrite=True)
                rv_mt = hl.read_matrix_table(rv_mt_path)

                if not rv_vcf_existence_outcome:
                    # remove fields not in the VCF
                    rv_mt = rv_mt.drop('gvcf_info')

                    # export to vcf rare variants only
                    export_vcf(rv_mt, rv_vcf_path)

            if not cv_vcf_existence_outcome or not cv_mt_existence_outcome:
                # common variants only
                cv_mt = mt.filter_rows(hl.min(mt.variant_qc.AF) >= cv_maf_threshold)

                if not cv_mt_existence_outcome:
                    # save chrom + rare variant mt for group file script
                    cv_mt.write(cv_mt_path, overwrite=True)

                cv_mt = hl.read_matrix_table(cv_mt_path)

                if not cv_vcf_existence_outcome:
                    # remove fields not in the VCF
                    cv_mt = cv_mt.drop('gvcf_info')

                    # export to vcf common variants only
                    export_vcf(cv_mt, cv_vcf_path)

        # check existence of index file (CV) separately
        cv_index_file_existence_outcome = can_reuse(f'{cv_vcf_path}.csi')
        logging.info(f'Does {cv_vcf_path}.csi exist? {cv_index_file_existence_outcome}')
        if not cv_index_file_existence_outcome:
            add_remove_chr_and_index_job(cv_vcf_path)

        # do the same for rare variants
        rv_index_file_existence_outcome = can_reuse(f'{rv_vcf_path}.csi')
        logging.info(f'Does {rv_vcf_path}.csi exist? {rv_index_file_existence_outcome}')
        if not rv_index_file_existence_outcome:
            add_remove_chr_and_index_job(rv_vcf_path)

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
