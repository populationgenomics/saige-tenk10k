#!/usr/bin/env python3
# pylint: disable=import-error,missing-module-docstring,no-value-for-parameter,wrong-import-position

__author__ = 'annacuomo'

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
import sys

import logging

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    # get_config,
    # init_batch,
    output_path,
)

import click

import pandas as pd

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
    expression_files_prefix: str,
    chromosome: int,
    cell_type: str,
    gene_info_tsv: str,
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
    # get single-cell expression for the cell type of interest
    expression_tsv_path = dataset_path(
        os.path.join(
            expression_files_prefix,
            'expression_files',
            f'{cell_type}_expression.tsv',
        )
    )
    expression_df = pd.read_csv(expression_tsv_path, sep='\t', index_col=0)
    # extract all genes
    all_genes = expression_df.columns.values
    # select only genes on relevant chromosome
    gene_info_df = pd.read_csv(gene_info_tsv, sep='\t')
    genes_chrom = gene_info_df[gene_info_df['chr'] == chromosome].index.values
    common_genes = set(all_genes).intersection(set(genes_chrom))
    # return expression for the correct chromosomes only
    return expression_df[:, common_genes]


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
    expression_df,
    cov_df,
    smf_df,
):
    """Combine files to build final input

    Input:
    - Expression dataframe
    - Covariates dataframe
    - Sample mapping file, mapping cells to donors
    """
    # cells_exprs_smf = set(expression_df.index.values).intersection(
    #     set(smf_df['cell'].unique())
    # )
    # common_cells = set(cov_df.index.values).intersection(cells_exprs_smf)
    pheno_cov_df = pd.concat([expression_df, cov_df, smf_df], axis=1)
    return pheno_cov_df


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--gene-info-tsv')
@click.option('--expression-files-prefix')
@click.option('--sample-mapping-file-path')
def expression_pipeline(
    gene_info_tsv: str,
    expression_files_prefix: str,
    celltypes: str,
    sample_mapping_file_path: str,
):
    """
    Run expression processing pipeline
    """
    celltype_list = celltypes.split(',')
    chromosome_list = celltypes.split(',')
    logging.info(f'Cell types to run: {celltype_list}')
    logging.info(f'Chromosomes to run: {chromosome_list}')

    smf_df = pd.read_csv(sample_mapping_file_path, sep='\t')

    for celltype in celltype_list:
        # get covariates (cell type specific)
        cov_df = get_celltype_covariates(
            expression_files_prefix=expression_files_prefix, cell_type=celltype
        )
        for chromosome in chromosome_list:
            # get expression (cell type + chromosome)
            expr_df = get_chrom_celltype_expression(
                chromosome=chromosome, cell_type=celltype, gene_info_tsv=gene_info_tsv
            )
            # combine files
            pheno_cov_df = build_pheno_cov_filename(
                cov_df=cov_df, expression_df=expr_df, smf_df=smf_df
            )

            # write to output
            pheno_cov_filename = to_path(
                output_path(f'input_files/{chromosome}_{celltype}.csv')
            )
            with pheno_cov_filename.open('w') as pcf:
                pheno_cov_df.to_csv(pcf, index=False)


if __name__ == '__main__':
    expression_pipeline()
