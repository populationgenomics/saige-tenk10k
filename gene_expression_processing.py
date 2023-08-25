#!/usr/bin/env python3
# pylint: disable=import-error

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

import pandas as pd

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