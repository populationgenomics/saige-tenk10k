#!/usr/bin/env python3

"""
This script will

- open anndata expression files
- create pheno_cov_files
- create cis window files

these files will be used as inputs for the
SAIGE-QTL association pipeline.
"""

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    get_batch,
    get_config,
    init_batch,
    output_path,
)
import click
import scanpy as sc

SCANPY_IMAGE = get_config()['images']['scanpy']


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--expression-files-prefix')
@click.option('--min-pct-expr', type=int, default=5)
@click.option('--cis-window-size', type=int, default=100000)
@click.option(
    '--max-gene-concurrency',
    type=int,
    default=50,
    help=(
        'To avoid resource starvation, set this concurrency to limit '
        'horizontal scale. Higher numbers have a better walltime, but '
        'risk jobs that are stuck (which are expensive)'
    ),
)
def main(
    celltypes: str,
    chromosomes: str,
    expression_files_prefix: str,
    min_pct_expr: int,
    cis_window_size: int,
    max_gene_concurrency=int,
):
    """
    Run expression processing pipeline
    """
    expression_h5ad_path = to_path(
        dataset_path(f'scrna-seq/CellRegMap_input_files/expression_objects/sce22.h5ad')
    ).copy('here.h5ad')
    expression_adata = sc.read(expression_h5ad_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
