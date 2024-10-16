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
   --dataset "bioheart" \
   --access-level "standard" \
   --output-dir saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/genotypes/ \
    python3 get_genotype_vcf.py --vds-path=gs://cpg-bioheart-main/vds/tenk10k1-0.vds --chromosomes chr2 \
    --relateds-to-drop-path=gs://cpg-bioheart-main-analysis/large_cohort/bioheart1-0/relateds_to_drop.ht
"""

import logging

import click
import pandas as pd

import hail as hl

from hail.methods import export_plink, export_vcf

from cpg_utils import to_path
from cpg_utils.config import get_config, output_path
from cpg_utils.hail_batch import get_batch, init_batch


VARS_PER_PARTITION = 20000


def can_reuse(path: str):
    if path and to_path(path).exists():
        return True
    return False


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


def remove_chr_from_bim(input_bim: str, output_bim: str, bim_renamed: str):
    """
    Method powered by Gemini

    Reads a PLINK .bim file, modifies the "chrom" column to numerical values,
    and saves the modified data to a new file.

    Args:
      input_bim: Path to the original .bim file.
      output_bim: Path to save the modified .bim file.
      bim_renamed: Path to success file.
    """
    # Read the .bim file into a DataFrame
    data = pd.read_csv(
        input_bim,
        sep='\t',
        header=None,
        names=['chrom', 'rsid', 'cm', 'bp', 'ref', 'alt'],
    )
    # Extract numerical chromosome values
    data['chrom'] = data['chrom'].str.extract('(\d+)')[0]

    # Save the modified DataFrame to a new .bim file
    with to_path(output_bim).open('w') as f:
        data.to_csv(f, sep='\t', header=None, index=False)

    # Write empty success file
    with to_path(bim_renamed).open('w') as fp:  # noqa
        pass


# inputs:
@click.option('--vds-path', help=' GCP gs:// path to the VDS')
@click.option('--chromosomes', help=' e.g., chr22,chrX ')
@click.option(
    '--relateds-to-drop-path',
    default='gs://cpg-bioheart-main-analysis/large_cohort/bioheart1-0/relateds_to_drop.ht',
)
@click.option('--cv-maf-threshold', default=0.01)
@click.option('--rv-maf-threshold', default=0.01)
@click.option('--vre-mac-threshold', default=20)
@click.option('--vre-n-markers', default=2000)
@click.option('--exclude-multiallelic', is_flag=False)
@click.option('--exclude-indels', is_flag=False)
@click.option('--plink-job-storage', default='1G')
@click.option('--ld-prune-r2', default=0.2)
@click.option('--ld-prune-bp-window-size', default=500000)
@click.command()
def main(
    vds_path: str,
    chromosomes: str,
    relateds_to_drop_path: str,
    cv_maf_threshold: float,
    rv_maf_threshold: float,
    vre_mac_threshold: int,
    vre_n_markers: int,
    exclude_multiallelic: bool,
    exclude_indels: bool,
    plink_job_storage: str,
    ld_prune_r2: float,
    ld_prune_bp_window_size: int,
):
    """
    Write genotypes as VCF
    """

    init_batch(worker_memory='highmem')

    vds = hl.vds.read_vds(vds_path)
    vds_name = vds_path.split('/')[-1].split('.')[0]

    for chromosome in chromosomes.split(','):
        # create paths and check if they exist already
        cv_vcf_path = output_path(
            f'vds-{vds_name}/{chromosome}_common_variants.vcf.bgz'
        )
        cv_vcf_existence_outcome = can_reuse(cv_vcf_path)
        logging.info(f'Does {cv_vcf_path} exist? {cv_vcf_existence_outcome}')

        rv_vcf_path = output_path(f'vds-{vds_name}/{chromosome}_rare_variants.vcf.bgz')
        rv_vcf_existence_outcome = can_reuse(rv_vcf_path)
        logging.info(f'Does {rv_vcf_path} exist? {rv_vcf_existence_outcome}')

        rv_mt_path = output_path(f'vds-{vds_name}/{chromosome}_rare_variants.mt')
        rv_mt_existence_outcome = can_reuse(rv_mt_path)
        logging.info(f'Does {rv_mt_path} exist? {rv_mt_existence_outcome}')

        if (
            not cv_vcf_existence_outcome
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
            # this will get dropped as the vds file will already be clean
            related_ht = hl.read_table(relateds_to_drop_path)
            related_samples = hl.literal(related_ht.s.collect())
            mt = mt.filter_cols(~related_samples.contains(mt['s']))

            # filter out loci & variant QC
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
                    rv_mt.write(rv_mt_path)
                    rv_mt = hl.read_matrix_table(rv_mt_path)

                if not rv_vcf_existence_outcome:
                    # remove fields not in the VCF
                    rv_mt = rv_mt.drop('gvcf_info')

                    # export to vcf rare variants only
                    export_vcf(rv_mt, rv_vcf_path)

            if not cv_vcf_existence_outcome:
                # common variants only
                cv_mt = mt.filter_rows(hl.min(mt.variant_qc.AF) >= cv_maf_threshold)

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

    # subset variants for variance ratio estimation
    vre_plink_path = output_path(f'vds-{vds_name}/vre_plink_2000_variants')
    vre_bim_path = f'{vre_plink_path}.bim'
    plink_existence_outcome = can_reuse(vre_bim_path)
    logging.info(f'Does {vre_bim_path} exist? {plink_existence_outcome}')
    if not plink_existence_outcome:
        # keep autosome chromosomes only
        vds = hl.vds.filter_chromosomes(vds, keep_autosomes=True)
        # split multiallelic loci pre densifying to mt
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)
        # densify to mt
        mt = hl.vds.to_dense_mt(vds)

        # drop a checkpoint here
        dense_checkpoint = output_path(
            'mt_from_dense_vds_checkpoint.mt', category='tmp'
        )

        if (to_path(dense_checkpoint) / '_SUCCESS').exists():
            print(f'Reading existing checkpoint from {dense_checkpoint}')
            mt = hl.read_matrix_table(dense_checkpoint)
        else:
            print(f'Writing new checkpoint to {dense_checkpoint}')
            mt = mt.checkpoint(dense_checkpoint)

        # filter out related samples from vre too
        # this will get dropped as the vds file will already be clean
        related_ht = hl.read_table(relateds_to_drop_path)
        related_samples = hl.literal(related_ht.s.collect())
        mt = mt.filter_cols(~related_samples.contains(mt['s']))

        # again filter for biallelic SNPs
        mt = mt.filter_rows(
            ~(mt.was_split)  # biallelic (exclude multiallelic)
            & (hl.len(mt.alleles) == 2)  # remove hom-ref
            & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
        )
        mt = hl.variant_qc(mt)

        # minor allele count (MAC) > {vre_mac_threshold}
        mt = mt.filter_rows(hl.min(mt.variant_qc.AC) > vre_mac_threshold)

        if (n_ac_vars := mt.count_rows()) == 0:
            raise ValueError('No variants left, exiting!')
        logging.info(f'MT filtered to common enough variants, {n_ac_vars} left')

        # drop a checkpoint here
        common_checkpoint = output_path('common_checkpoint.mt', category='tmp')

        if (to_path(common_checkpoint) / '_SUCCESS').exists():
            print(f'Reading existing checkpoint from {common_checkpoint}')
            mt = hl.read_matrix_table(common_checkpoint)
        else:
            print(f'Writing new checkpoint to {common_checkpoint}')
            mt = mt.checkpoint(common_checkpoint)

        logging.info(f'common checkpoint written, MT size: {mt.count()}')
        mt.describe()

        # since pruning is very costly, subset first a bit
        if n_ac_vars > (vre_n_markers * 100):
            logging.info('subset completed')
            mt = mt.sample_rows(p=0.01, seed=0)

        # set a checkpoint, and either re-use or write
        post_downsampling_checkpoint = output_path(
            'common_subset_checkpoint.mt', category='tmp'
        )

        # overwrite=True to force re-write, requires full permissions
        mt = mt.checkpoint(post_downsampling_checkpoint, overwrite=True)

        # perform LD pruning
        pruned_variant_table = hl.ld_prune(
            mt.GT, r2=ld_prune_r2, bp_window_size=ld_prune_bp_window_size
        )

        # write the resulting table
        pruned_variant_table = pruned_variant_table.checkpoint(
            output_path('pruned_result.ht', category='tmp')
        )

        # filter MT to only include pruned variants
        mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

        # only count once
        remaining_variants = mt.count_rows()

        logging.info(f'pruning completed, {remaining_variants} variants left')
        # randomly sample {vre_n_markers} variants
        mt = mt.sample_rows((vre_n_markers * 1.1) / remaining_variants, seed=0)
        mt = mt.head(vre_n_markers)

        # export to plink common variants only for sparse GRM
        export_plink(mt, vre_plink_path, ind_id=mt.s)
        logging.info('plink export completed')

    # success file for chr renaming in bim file
    bim_renamed_path = output_path(f'vds-{vds_name}/bim_renamed.txt')
    bim_renamed_existence_outcome = can_reuse(bim_renamed_path)
    logging.info(
        f'Have the chr been renamed in the bim file? {bim_renamed_existence_outcome}'
    )
    if not bim_renamed_existence_outcome:
        # saige requires numerical values for chromosomes, so
        # removing "chr" from the bim file
        plink_input_bim = get_batch().read_input(vre_bim_path)
        remove_chr_job = get_batch().new_python_job(name='remove chr from plink bim')
        remove_chr_job.cpu(1)
        remove_chr_job.storage(plink_job_storage)
        # remove chr, then write direct to the BIM source location
        remove_chr_job.call(
            remove_chr_from_bim, plink_input_bim, vre_bim_path, bim_renamed_path
        )
        logging.info('chr removed from bim')

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
