#!/usr/bin/env python3

"""
This script will attempt to extract genotype PCs
for TOB-WGS + BioHEART


To run:

analysis-runner \
    --description "extract geno PCs" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/covariates/vds-input" \
    python3 get_geno_pcs.py --vds-version 1-0

"""

import click
import hail as hl

from cpg_utils.hail_batch import dataset_path, init_batch, output_path


@click.option('--vds-version', help=' e.g., 1-0 ')
@click.command()
def main(vds_version):
    """
    Obtain genotype PCs and write loading and scores
    """

    init_batch(worker_memory='highmem')

    # load vds object
    vds_path = dataset_path(f'vds/{vds_version}.vds')
    vds = hl.vds.read_vds(vds_path)
        # get pcs
    _, scores, loadings = hl.hwe_normalized_pca(vds.LGT, k=15, compute_loadings=True)
    # write
    scores.write(output_path(f'{vds_version}/geno_pca_scores.ht'), overwrite=True)
    loadings.write(output_path(f'{vds_version}/geno_pca_loadings.ht'), overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
