#!/usr/bin/env python3

"""
This script will

- extract common variants
- export as VCF (one per chrom)
- also create a plink file for a subset of variants
for variance ratio estimation

this will be used as input for the
SAIGE-QTL association pipeline.


To run:

In test:

analysis-runner \
    --description "get common variant VCF" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/genotypes/" \
    python3 get_genotype_vcf.py --vds-version 1-0 --chromosomes chr1,chr2,chr22 --vre-mac-threshold 1

In main:

analysis-runner \
    --description "get common variant VCF" \
    --dataset "bioheart" \
    --access-level "full" \
    --output-dir "saige-qtl/input_files/genotypes/" \
    python3 get_genotype_vcf.py --vds-version v1-0 --chromosomes chr1,chr2,chr22

"""

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, init_batch, output_path

import click
import logging
import random
import pandas as pd

import hail as hl
from hail.methods import export_plink, export_vcf


# this is a guess, let's see how it performs...
VARS_PER_PARTITION = 20000


def checkpoint_mt(mt: hl.MatrixTable, checkpoint_path: str, force: bool = False):
    """
    Checkpoint a MatrixTable to a file.
    If the checkpoint exists, read instead

    inspired by this thread
    https://discuss.hail.is/t/best-way-to-repartition-heavily-filtered-matrix-tables/2140

    Args:
      mt: MatrixTable to checkpoint.
      checkpoint_path: Path to save the MatrixTable to.
      force: Whether to overwrite an existing checkpoint
    """

    # create a temp checkpoint path
    temp_checkpoint_path = checkpoint_path + '.temp'
    logging.info(f'Checkpoint temp: {temp_checkpoint_path}')

    # either force an overwrite, or write the first version
    if force or not to_path(temp_checkpoint_path).exists():
        logging.info(f'Writing new temp checkpoint to {temp_checkpoint_path}')
        mt = mt.checkpoint(temp_checkpoint_path, overwrite=True)

    elif (
        to_path(checkpoint_path).exists()
        and (to_path(checkpoint_path) / '_SUCCESS').exists()
    ):
        logging.info(f'Reading non-temp checkpoint from {checkpoint_path}')
        return hl.read_matrix_table(checkpoint_path)

    # unless forced, if the data exists, read it
    elif (
        to_path(temp_checkpoint_path).exists()
        and (to_path(temp_checkpoint_path) / '_SUCCESS').exists()
    ):
        logging.info(f'Reading existing temp checkpoint from {temp_checkpoint_path}')
        mt = hl.read_matrix_table(temp_checkpoint_path)

    else:
        raise FileNotFoundError('Checkpoint exists but is incomplete, was not forced')

    logging.info(
        f'Dimensions of MT: {mt.count()}, across {mt.n_partitions()} partitions'
    )

    # repartition to a reasonable number of partitions, then re-write
    hl.read_matrix_table(
        temp_checkpoint_path, _n_partitions=mt.count_rows() // VARS_PER_PARTITION or 1
    ).write(checkpoint_path)

    # delete the temp checkpoint
    hl.current_backend().fs.rmtree(temp_checkpoint_path)

    return mt


def can_reuse(path: str):
    if path and to_path(path).exists():
        return True
    return False


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
@click.option('--vds-version', help=' e.g., 1-0 ')
@click.option('--chromosomes', help=' e.g., chr22,chrX ')
@click.option('--cv-maf-threshold', default=0.01)
@click.option('--vre-mac-threshold', default=20)
@click.option('--vre-n-markers', default=2000)
@click.option('--exclude-multiallelic', is_flag=False)
@click.option('--exclude-indels', is_flag=False)
@click.option('--plink-job-storage', default='1G')
@click.command()
def main(
    vds_version: str,
    chromosomes: str,
    cv_maf_threshold: float,
    vre_mac_threshold: int,
    vre_n_markers: int,
    exclude_multiallelic: bool,
    exclude_indels: bool,
    plink_job_storage: str,
):
    """
    Write genotypes as VCF
    """

    init_batch(worker_memory='highmem')

    vds_path = dataset_path(f'vds/{vds_version}.vds')
    vds = hl.vds.read_vds(vds_path)

    for chromosome in chromosomes.split(','):

        # create path and check if it exists already
        cv_vcf_path = output_path(
            f'vds{vds_version}/{chromosome}_common_variants.vcf.bgz'
        )
        vcf_existence_outcome = can_reuse(cv_vcf_path)
        logging.info(f'Does {cv_vcf_path} exist? {vcf_existence_outcome}')
        if not vcf_existence_outcome:

            # consider only relevant chromosome
            chrom_vds = hl.vds.filter_chromosomes(vds, keep=chromosome)

            # split multiallelic loci (necessary pre-densifying)
            chrom_vds = hl.vds.split_multi(chrom_vds, filter_changed_loci=True)

            # densify to matrix table object
            mt = hl.vds.to_dense_mt(chrom_vds)

            # filter out loci & variant QC
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # remove hom-ref
            if exclude_multiallelic:  # biallelic only (exclude multiallelic)
                mt = mt.filter_rows(~(mt.was_split))
            if exclude_indels:  # SNPs only (exclude indels)
                mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
            # logging.info(f'Number of variants left after filtering: {mt.count()}')

            mt = hl.variant_qc(mt)

            # common variants only
            cv_mt = mt.filter_rows(mt.variant_qc.AF[1] > cv_maf_threshold)
            # logging.info(f'Number of common variants left: {cv_mt.count()}')

            # remove fields not in the VCF
            cv_mt = cv_mt.drop('gvcf_info')

            # export to vcf common variants only
            export_vcf(cv_mt, cv_vcf_path)

        # check existence of index file separately
        index_file_existence_outcome = can_reuse(f'{cv_vcf_path}.csi')
        logging.info(f'Does {cv_vcf_path}.csi exist? {index_file_existence_outcome}')
        if not index_file_existence_outcome:
            # remove chr & add index file using bcftools
            vcf_input = get_batch().read_input(cv_vcf_path)

            vcf_size = to_path(cv_vcf_path).stat().st_size
            storage_required = ((vcf_size // 1024**3) or 1) * 2.2
            bcftools_job = get_batch().new_job(name='remove chr and index vcf')
            bcftools_job.declare_resource_group(
                output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.csi': '{root}.vcf.bgz.csi'}
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
                bcftools annotate --rename-chrs chr_update.txt {vcf_input} | \\
                bgzip -c > {bcftools_job.output['vcf.bgz']}
                bcftools index -c {bcftools_job.output['vcf.bgz']}
            """
            )
            logging.info('VCF rename/index jobs scheduled!')

            # save both output files
            get_batch().write_output(bcftools_job.output, cv_vcf_path.removesuffix('.vcf.bgz'))

    # subset variants for variance ratio estimation
    vre_plink_path = output_path(f'vds{vds_version}/vre_plink_2000_variants')
    vre_bim_path = f'{vre_plink_path}.bim'
    plink_existence_outcome = can_reuse(vre_bim_path)
    logging.info(f'Does {vre_bim_path} exist? {plink_existence_outcome}')
    if not plink_existence_outcome:
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)
        mt = hl.vds.to_dense_mt(vds)

        # set a checkpoint, and either re-use or write
        post_dense_checkpoint = output_path('post_dense_checkpoint.mt', category='tmp')
        mt = checkpoint_mt(mt, post_dense_checkpoint)

        # again filter for biallelic SNPs
        mt = mt.filter_rows(
            ~(mt.was_split)  # biallelic (exclude multiallelic)
            & (hl.len(mt.alleles) == 2)  # remove hom-ref
            & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
        )
        mt = hl.variant_qc(mt)

        # minor allele count (MAC) > {vre_mac_threshold}
        vre_mt = mt.filter_rows(mt.variant_qc.AC[1] > vre_mac_threshold)

        # set a checkpoint, and either re-use or write
        post_common_checkpoint = output_path(
            'common_reduced_checkpoint.mt', category='tmp'
        )
        vre_mt = checkpoint_mt(vre_mt, post_common_checkpoint)

        if (n_ac_vars := vre_mt.count_rows()) == 0:
            raise ValueError('No variants left, exiting!')
        logging.info(f'MT filtered to common enough variants, {n_ac_vars} left')

        # since pruning is very costly, subset first a bit
        random.seed(0)
        vre_mt = vre_mt.sample_rows(p=0.01)
        logging.info('subset completed')

        # perform LD pruning
        pruned_variant_table = hl.ld_prune(vre_mt.GT, r2=0.2, bp_window_size=500000)
        vre_mt = vre_mt.filter_rows(hl.is_defined(pruned_variant_table[vre_mt.row_key]))

        post_prune_checkpoint = output_path('post_prune_checkpoint.mt', category='tmp')
        vre_mt = checkpoint_mt(vre_mt, post_prune_checkpoint)

        logging.info(f'pruning completed, {vre_mt.count_rows()} variants left')
        # randomly sample {vre_n_markers} variants
        random.seed(0)
        vre_mt = vre_mt.sample_rows((vre_n_markers * 1.1) / vre_mt.count_rows())
        vre_mt = vre_mt.head(vre_n_markers)
        logging.info(f'sampling completed, {vre_mt.count()} variants left')

        # export to plink common variants only for sparse GRM
        export_plink(vre_mt, vre_plink_path, ind_id=vre_mt.s)
        logging.info('plink export completed')

    # success file for chr renaming in bim file
    bim_renamed_path = output_path(f'vds{vds_version}/bim_renamed.txt')
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
