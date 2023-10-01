#!/usr/bin/env python3


"""
Hail Batch workflow to create gene expression files.
This script will:

- select genes to test based on expression (per cell type)
- build chromosome & cell type specific phenotype covariate files
- use gene info to create cis-window files
- export pheno_cov files as tsv files
- export cis_window files as tsv files

More details in README
output files in tob_wgs_genetics/saige_qtl/input

analysis-runner \
    --dataset tob-wgs \
    --access-level test \
    --output-dir 'tob_wgs_genetics/saige_qtl/hope-test-input' \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:0.0.3 \
    --description 'scRNA-seq processing batch job test' \
    python3 gene_expression_processing_simple.py \
    --celltypes=B_IN --chromosomes=chr22 \
    --gene-info-tsv=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv \
    --sample-mapping-file-path=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv \
    --expression-files-prefix=hope-test

"""

import math

import click
import hail as hl
import hailtop.batch as hb
import pandas as pd
import scanpy as sc

from cpg_utils import to_path

from cpg_utils.hail_batch import copy_common_env, dataset_path, output_path
from cpg_workflows.batch import get_batch



def get_chrom_celltype_expression(
    gene_info_df,
    expression_files_prefix: str,  # tob_wgs_genetics/saige_qtl/input/
    chromosome: str,
    cell_type: str,
) -> sc.AnnData:
    """Extracts relevant expression info

    Input:
    - chromosome & cell type of interest
    - path to (single-cell) expression files, one tsv file per cell type,
    rows are cells and columns are genes
    - path to dataframe containing gene info, for each gene (=row),
    specifies chrom, start, end and strand

    Output: expression adata object for only relevant genes
    """

    # first line is where the file is now,
    # but (second line) the files will eventually be in the below folder
    # and split by cell type (at least this all naive B cells only)
    expression_h5ad_path = to_path(
        dataset_path(
            f'scrna-seq/CellRegMap_input_files/expression_objects/sce22.h5ad'
        )
    ).copy('here.h5ad')
    
    expression_adata = sc.read(expression_h5ad_path)

    # select only genes on relevant chromosome
    genes_chrom = gene_info_df[gene_info_df['chr'] == chromosome].gene_name
    # return expression for the correct chromosomes only
    return expression_adata[:, expression_adata.var_names.isin(genes_chrom)]

@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--gene-info-tsv')
@click.option('--expression-files-prefix')
@click.option('--sample-mapping-file-path')
@click.option('--min-pct-expr', type=int, default =5)
@click.option('--cis-window-size', type=int,default=100000)
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
def expression_pipeline(
    celltypes: str,
    chromosomes: str,
    gene_info_tsv: str,
    expression_files_prefix: str,
    sample_mapping_file_path: str,
):
    """
    Run expression processing pipeline
    """

    b = get_batch()
    # create phenotype covariate files
    smf_df = pd.read_csv(sample_mapping_file_path, sep='\t')
    gene_info_df = pd.read_csv(gene_info_tsv, sep='\t')

    
    # set jobs running
    b.run(wait=False)


if __name__ == '__main__':
    expression_pipeline()  # pylint: disable=no-value-for-parameter
