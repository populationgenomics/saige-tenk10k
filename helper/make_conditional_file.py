#!/usr/bin/env python3
"""
This script summarise SAIGE-QTL single-variant test
raw p-values from common (or rare) variant results
and plots QQ plots and histograms of the p-values

To run:

analysis-runner \
    --description "make conditional files" \
    --dataset "bioheart" \
    --access-level "test" \
    --image 'australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16' \
    --output-dir "saige-qtl/" \
    python3 plotter/summarise_and_qq_plotter.py \
        --celltype='B_naive' \
        --results-path=gs://cpg-bioheart-main-analysis/saige-qtl/bioheart_n990_and_tob_n1055/output_files/sample_perm0/output_files/ \
        --title='shuffled'
"""

import click
import matplotlib.pyplot as plt
import hail as hl

import numpy as np
import pandas as pd
from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch, output_path

from multipy.fdr import qvalue


@click.command()
@click.option('--celltypes', required=True, help='separated by comma')
@click.option('--qtl-results-path', required=True)
@click.option('--conditional-files-output-path', required=True)
@click.option('--add-to-existing-file', default='nope')
@click.option('--qv-significance-threshold', default=0.05)
def make_condition_file(
    celltypes: str,
    results_path: str,
    conditional_files_output_path: str,
    add_to_existing_file: str,
    qv_significance_threshold: float,
):
    """
    gets summarised results from a run of SAIGE-QTL CV (step3)
    extracts significant genes and corresponding top variants
    builds conditional files (per celltype)
    """
    init_batch()

    for celltype in celltypes.split(','):

        gene_level_results = (
            f'{results_path}/summary_stats/{celltype}_all_cis_cv_gene_level_results.tsv'
        )
        gene_level_df = pd.read_csv(gene_level_results, sep='\t')
        _, qvals = qvalue(gene_level_df['ACAT_p'])
        gene_level_df['qvalue'] = qvals
        significant_gene_df = gene_level_df[
            gene_level_df['qvalue'] < qv_significance_threshold
        ]
        if add_to_existing_file == 'nope':
            conditional_df = significant_gene_df
            conditional_df['variants_to_condition_on'] = conditional_df['top_MarkerID']
        elif add_to_existing_file != 'nope':
            old_conditional_df = pd.read_csv(add_to_existing_file, sep='\t')
            # open file (called 'add_to_existing_file')
            # subset to genes in current significant_gene_df
            # add top marker from current significant_gene_df to the condition in those
            # reorder by genomic coordinates
        # write conditional_df to file
        conditional_file = (
            f'{conditional_files_output_path}/{celltype}_conditional_file.tsv'
        )
        conditional_df.to_csv(conditional_file, sep='\t')


if __name__ == '__main__':
    make_condition_file()
