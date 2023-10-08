#!/usr/bin/env python3
# pylint: disable=broad-exception-raised,import-error,import-outside-toplevel,missing-module-docstring,no-value-for-parameter,too-many-arguments,too-many-branches,too-many-locals,too-many-statements,wrong-import-order,wrong-import-position


"""
This script will:

- export genotype subset as plink files for fitting to a null model

More details in README
output files in tob_wgs_genetics/saige_qtl/input
"""

import logging
import random
from cpg_utils.hail_batch import (
    get_config,
    init_batch,
    output_path,
)

import click
import hail as hl

HAIL_IMAGE = get_config()['images']['hail']

# only needs to be run once for a given cohort (e.g., OneK1K / TOB)
def ld_prune(
    mt_path: str,  # 'tob_wgs/filtered.mt'
    vre_plink_path: str,  # 'tob_wgs/vr_plink_2000_variants
    vre_mac_threshold: int = 20,
    vre_n_markers: int = 2000,
):
    """Subset hail matrix table

    Outputs: 
    plink file containing a random subset of 2,000 variants that are not in LD 
    """
    from hail.methods import export_plink

    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # subset variants for variance ratio estimation
    # minor allele count (MAC) > {vre_n_markers}
    vre_mt = mt.filter_rows(mt.variant_qc.AC[0] > vre_mac_threshold)
    n_ac_vars = vre_mt.count()[0]  # to avoid evaluating this 2X
    logging.info(f'Number of variants post AC filter: {n_ac_vars}')
    if n_ac_vars == 0:
        logging.info('No variants left, exit')
        return
    # since pruning is very costly, subset first a bit
    random.seed(0)
    vre_mt = vre_mt.sample_rows(p=0.01)
    logging.info(f'Initial subset of variants: {vre_mt.count()[0]}')
    # perform LD pruning
    pruned_variant_table = hl.ld_prune(vre_mt.GT, r2=0.2, bp_window_size=500000)
    vre_mt = vre_mt.filter_rows(hl.is_defined(pruned_variant_table[vre_mt.row_key]))
    logging.info(f'Subset of variants after pruning: {vre_mt.count()[0]}')
    # randomly sample {vre_n_markers} variants
    random.seed(0)
    vre_mt = vre_mt.sample_rows((vre_n_markers * 1.1) / vre_mt.count()[0])
    vre_mt = vre_mt.head(vre_n_markers)
    logging.info(f'Subset to {vre_n_markers} variants: {vre_mt.count()[0]}')

    # export to plink common variants only for sparse GRM
    export_plink(vre_mt, vre_plink_path, ind_id=vre_mt.s)


@click.command()
@click.option('--mt-path')

def main(
    mt_path: str
):
    
    vre_plink_path = output_path('vr_plink_20k_variants')

    ld_prune(
        mt_path=mt_path,
        vre_plink_path=vre_plink_path,
    )


if __name__ == '__main__':
    main()
