#!/usr/bin/env python3

"""
This script will extract a subset of n independent common variants
for variance ratio estimation, where
- n is user specified, 2,000 by default
- independent is intended as after pruning, with a user specified LD threshold (0.2 by default)
- common is user specified, by default MAC > 20

this will be used as input for the
SAIGE-QTL association pipeline.

analysis-runner \
   --description "get variance ratio estimation plink files" \
   --dataset "tenk10k" \
   --access-level "standard" \
   --output-dir saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/genotypes/ \
    python3 get_vre_plink.py --vds-path=gs://cpg-tenk10k-main/vds/tenk10k-genome-2-0.vds \
    --relateds-to-drop-path=gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/relateds_to_drop.ht \
    --qc-pass-variants-path=gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/variants_qc.ht
"""

import logging

import click
import pandas as pd

import hail as hl

from hail.methods import export_plink

from cpg_utils import to_path
from cpg_utils.config import output_path
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
@click.option(
    '--relateds-to-drop-path',
    default='gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/relateds_to_drop.ht',
)
@click.option(
    '--qc-pass-variants-path',
    default='gs://cpg-tenk10k-main/large_cohort/tenk10k-genome-2-3-eur/variants_qc.ht',
)
@click.option('--vre-mac-threshold', default=20)
@click.option('--vre-n-markers', default=2000)
@click.option('--plink-job-storage', default='1G')
@click.option('--ld-prune-r2', default=0.2)
@click.option('--ld-prune-bp-window-size', default=500000)
@click.command()
def main(
    vds_path: str,
    relateds_to_drop_path: str,
    qc_pass_variants_path: str,
    vre_mac_threshold: int,
    vre_n_markers: int,
    plink_job_storage: str,
    ld_prune_r2: float,
    ld_prune_bp_window_size: int,
):
    """
    Write one subset of the genotypes
    as plink files
    """

    init_batch(worker_memory='highmem')

    vds = hl.vds.read_vds(vds_path)
    vds_name = vds_path.split('/')[-1].split('.')[0]

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
            # force overwrite, in case a partially written checkpoint existed
            mt = mt.checkpoint(dense_checkpoint, overwrite=True)

        # filter out related samples from vre too
        related_ht = hl.read_table(relateds_to_drop_path)
        related_samples = hl.literal(related_ht.s.collect())
        mt = mt.filter_cols(~related_samples.contains(mt['s']))

        # again filter for QC-pass, biallelic SNPs
        qc_pass_variants = hl.read_table(qc_pass_variants_path)
        mt = mt.semi_join_rows(qc_pass_variants)
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
            mt = mt.checkpoint(common_checkpoint, overwrite=True)

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
