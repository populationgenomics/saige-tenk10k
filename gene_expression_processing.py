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
"""

# import python modules
import os
import logging
import sys

import click
import hail as hl
import pandas as pd
import scanpy as sc

from cloudpathlib import AnyPath

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    output_path,
)


# use logging to print statements, display at info level
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr,
)


# adapted from https://github.com/populationgenomics/tob-wgs/blob/get-variants/
# scripts/eqtl_hail_batch/generate_eqtl_spearman.py#L34-L60
def filter_lowly_expressed_genes(expression_df, min_pct=5):
    """Remove genes with low expression across cells

    Input:
    expression_df: a data frame with cells as rows and genes as columns,
    containing normalised expression values (i.e., the number of molecules
    for each gene detected in each cell, normalised by sequencing depth).

    Returns:
    A filtered version of the input data frame, after removing columns (genes)
    expressed in less than {min_pct}% of the rows (cells).
    """

    # Remove genes with 0 expression in all cells
    expression_df = expression_df.loc[:, (expression_df != 0).any(axis=0)]
    genes_not_equal_zero = expression_df.iloc[:, 1:].values != 0
    n_expr_over_zero = pd.DataFrame(genes_not_equal_zero.sum(axis=0))
    percent_expr_over_zero = (n_expr_over_zero / len(expression_df.index)) * 100
    percent_expr_over_zero.index = expression_df.columns[1:]

    # Filter genes with less than {min_pct} percent cells with non-zero expression
    atleastminpercent = percent_expr_over_zero[(percent_expr_over_zero > min_pct)[0]]
    sample_ids = expression_df['sampleid']
    expression_df = expression_df[atleastminpercent.index]
    expression_df.insert(loc=0, column='sampleid', value=sample_ids)

    return expression_df


def get_chrom_celltype_expression(
    gene_info_df,
    # expression_files_prefix: str,
    chromosome: str,
    # cell_type: str,
):
    """Extracts relevant expression info

    Input:
    - chromosome & cell type of interest
    - path to (single-cell) expression files, one tsv file per cell type,
    rows are cells and columns are genes
    - path to dataframe containing gene info, for each gene (=row),
    specifies chrom, start, end and strand

    Output: expression dataframe for only relevant genes
    """
    # get single-cell expression for the cell type
    # and chromosome of interest (check)
    # expression_tsv_path = dataset_path(
    #     os.path.join(
    #         expression_files_prefix,
    #         'expression_files',
    #         cell_type,
    #         f'{chromosome}_expression.tsv',
    #     )
    # )
    # expression_df = pd.read_csv(expression_tsv_path, sep='\t', index_col=0)

    # this is where the file is now, but the files will eventually be in a different folder
    # and split by cell type (at least this all naive B cells only)
    expression_h5ad_path = AnyPath(
        dataset_path(
            f'scrna-seq/CellRegMap_input_files/expression_objects/sce{chromosome}.h5ad'
        )
    ).copy('here.h5ad')
    expression_adata = sc.read(expression_h5ad_path)

    # extract all genes
    all_genes = expression_adata.raw.var.index
    # all_genes = expression_df.columns.values
    # select only genes on relevant chromosome
    genes_chrom = gene_info_df[gene_info_df['chr'] == chromosome].index.values
    common_genes = set(all_genes).intersection(set(genes_chrom))
    # return expression for the correct chromosomes only
    # return expression_df[:, common_genes]
    return expression_adata[:, common_genes]  # check syntax


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
    covs_tsv_path = dataset_path(
        os.path.join(
            expression_files_prefix,
            'covariate_files',
            f'{cell_type}_covs.tsv',
        )
    )
    covs_df = pd.read_csv(covs_tsv_path, sep='\t', index_col=0)
    return covs_df


def build_pheno_cov_filename(
    gene_name,
    expression_adata,
    cov_df,
    smf_df,  # sample mapping file
):
    """Combine files to build final input
    reduce to one file per gene

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
    return pheno_cov_df


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
    gene_cis_df = pd.DataFrame(data, index=gene)
    return gene_cis_df


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--gene-info-tsv')
@click.option('--expression-files-prefix')
@click.option('--sample-mapping-file-path')
@click.option('--min-pct-expr')
@click.option('--cis-window-size')
def expression_pipeline(
    celltypes: str,
    chromosomes: str,
    gene_info_tsv: str,
    expression_files_prefix: str,
    sample_mapping_file_path: str,
    min_pct_expr: int = 5,
    cis_window_size: int = 100000,
):
    """
    Run expression processing pipeline
    """
    celltype_list = celltypes.split(',')
    chromosome_list = chromosomes.split(',')
    logging.info(f'Cell types to run: {celltype_list}')
    logging.info(f'Chromosomes to run: {chromosome_list}')

    # create phenotype covariate files
    smf_df = pd.read_csv(sample_mapping_file_path, sep='\t')
    gene_info_df = pd.read_csv(gene_info_tsv, sep='\t')

    for celltype in celltype_list:
        # get covariates (cell type specific)
        cov_df = get_celltype_covariates(
            expression_files_prefix=expression_files_prefix, cell_type=celltype
        )
        for chromosome in chromosome_list:
            # get expression (cell type + chromosome)
            expr_df = get_chrom_celltype_expression(
                gene_info_df=gene_info_df,
                expression_files_prefix=expression_files_prefix,
                chromosome=chromosome,
                cell_type=celltype,
            )
            # remove lowly expressed genes
            expr_df = filter_lowly_expressed_genes(
                expression_df=expr_df, min_pct=min_pct_expr
            )
            # combine files
            pheno_cov_df = build_pheno_cov_filename(
                cov_df=cov_df, expression_df=expr_df, smf_df=smf_df
            )

            # write to output
            pheno_cov_filename = to_path(
                output_path(f'input_files/pheno_cov_files/{chromosome}_{celltype}.csv')
            )
            with pheno_cov_filename.open('w') as pcf:
                pheno_cov_df.to_csv(pcf, index=False)

    # create gene cis window files
    for gene in gene_info_df.index.values:
        gene_cis_filename = to_path(
            output_path(f'input_files/cis_window_files/{gene}_{cis_window_size}bp.csv')
        )
        gene_cis_df = get_gene_cis_file(
            gene_info_df=gene_info_df,
            gene=gene,
            window_size=cis_window_size,
        )
        with gene_cis_filename.open('w') as gcf:
            gene_cis_df.to_csv(gcf, index=False)


if __name__ == '__main__':
    expression_pipeline()  # pylint: disable=no-value-for-parameter
