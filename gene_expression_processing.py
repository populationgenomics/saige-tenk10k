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

   

    analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files" \
    --access-level "test" \
    --output-dir "tob_wgs_genetics/saige_qtl/hope-test-input" \
     gene_expression_processing.py  --celltypes=B_IN --chromosomes=chr22 \
    --gene-info-tsv=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv \
    --sample-mapping-file-path=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv \
    --expression-files-prefix=hope-test

"""

import os
import logging
import math
import sys
#from tkinter.tix import CELL

import click
import hail as hl
import pandas as pd
from cpg_workflows.batch import get_batch

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

from cpg_utils.hail_batch import copy_common_env, dataset_path, output_path

config = get_config()


CELLREGMAP_IMAGE = (
    'australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:0.0.3'
)


def filter_lowly_expressed_genes(expression_adata, min_pct=5):
    import scanpy 
    """Remove genes with low expression across cells

    Input: adata with all genes

    Output: adata filtered
    """
    n_all_cells = len(expression_adata.obs.index)
    min_cells = math.ceil((n_all_cells * min_pct) / 100)
    expression_adata = scanpy.pp.filter_genes(expression_adata, min_cells=min_cells)
    assert isinstance(expression_adata, scanpy.AnnData)

    return expression_adata


def get_chrom_celltype_expression(
    gene_info_df,
    expression_files_prefix: str,  # tob_wgs_genetics/saige_qtl/input/
    chromosome: str,
    cell_type: str,
):
    import scanpy
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
    expression_h5ad_path = 'gs://cpg-tob-wgs-test/scrna-seq/CellRegMap_input_files/expression_objects/sce22.h5ad'
    expression_adata = scanpy.read(expression_h5ad_path)

    # select only genes on relevant chromosome
    genes_chrom = gene_info_df[gene_info_df['chr'] == chromosome].gene_name
    # return expression for the correct chromosomes only
    return expression_adata[:, expression_adata.var_names.isin(genes_chrom)]


def get_celltype_covariates(
    expression_files_prefix: str,
    cell_type: str,
):
    """Obtain cell type specific covariates

    Input:
    - cell type of interest
    - covariate files prefix

    Output: covariate df for cell type of interest
    """
    covs_tsv_path = dataset_path('tob_wgs_genetics/saige_qtl/input/covariate_chr22_B_IN_tester.csv')
    covs_df = pd.read_csv(covs_tsv_path, sep=',', index_col=0)
    return covs_df


def build_pheno_cov_filename(
    gene_name, expression_adata, cov_df, smf_df, out_path: str  # sample mapping file
):
    """
    Combine files to build final input
    reduce to one file per gene
    write the output to disc within this method

    Input:
    - Expression dataframe
    - Covariates dataframe
    - Sample mapping file, mapping cells to donors
    """
    gene_adata = expression_adata[:, gene_name]
    gene_mat = gene_adata.raw.X.todense()
    expression_df = pd.DataFrame(
        data=gene_mat.T, index=gene_adata.raw.var.index, columns=gene_adata.obs.index
    )
    pheno_cov_df = pd.concat([expression_df, cov_df, smf_df], axis=1)

    with to_path(out_path).open('w') as pcf:
        pheno_cov_df.to_csv(pcf, index=False)


def get_gene_cis_file(gene_info_df, gene: str, window_size: int):
    """Get gene cis window file"""
    # select the gene from df
    gene_info_gene = gene_info_df[gene_info_df['gene'] == gene]
    # get gene chromosome
    chrom = gene_info_gene['chr']
    # get gene body position (start and end) and add window
    left_boundary = max(1, int(gene_info_gene['start']) - window_size)
    right_boundary = min(
        int(gene_info_gene['end']) + window_size,
        hl.get_reference('GRCh38').lengths[chrom],
    )
    data = {'chromosome': chrom, 'start': left_boundary, 'end': right_boundary}
    # check if I need an index at all
    return pd.DataFrame(data, index=[gene])


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
    min_pct_expr: int,
    cis_window_size: int,
    max_gene_concurrency=int
):
    """
    Run expression processing pipeline
    """
    
    b =get_batch()

    logging.info(f'Cell types to run: {celltypes}')
    logging.info(f'Chromosomes to run: {chromosomes}')

    # create phenotype covariate files
    smf_df = pd.read_csv(sample_mapping_file_path, sep='\t')
    gene_info_df = pd.read_csv(gene_info_tsv, sep='\t')

    for celltype in celltypes.split(','):
        # get covariates (cell type specific)
        cov_df = get_celltype_covariates(
            expression_files_prefix=expression_files_prefix, cell_type=celltype
        )
        for chromosome in chromosomes.split(','):
            # get expression (cell type + chromosome)
            j = b.new_python_job(name='Get expression (cell type + chr)')
            j.image(CELLREGMAP_IMAGE)
            expr_adata = j.call(get_chrom_celltype_expression,gene_info_df,expression_files_prefix,chromosome,celltype)
            #f = b.new_python_job(name = 'remove lowly expressed genes')
           # filter_adata = f.call(filter_lowly_expressed_genes,expression_adata=expr_adata, min_pct=min_pct_expr)

    b.run(wait=False)

if __name__ == '__main__':
    expression_pipeline()  # pylint: disable=no-value-for-parameter
