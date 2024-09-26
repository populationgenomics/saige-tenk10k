#!/usr/bin/env python3
"""
This script summarise SAIGE-QTL single-variant test
raw p-values from common (or rare) variant results
and plots QQ plots and histograms of the p-values

To run:

analysis-runner \
    --description "plot p-values" \
    --dataset "bioheart" \
    --access-level "test" \
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


@click.command()
@click.option('--celltypes', required=True, help='separated by comma')
@click.option('--results-path', required=True)
@click.option('--title', default='shuffled')
@click.option('--common-or-rare', default='common')
def plot_pvalues(
    celltypes: str,
    results_path: str,
    title: str,
    common_or_rare: str,
):
    """
    combines the results for a given cell type,
    saving both all results and top SNP per gene results
    plots both a histogram and a QQ plot of the association p-values
    """
    init_batch()

    for celltype in celltypes.split(','):

        # collect all raw p-value files
        if common_or_rare == 'common':
            existing_assoc_results = [
                str(file)
                for file in to_path(results_path).glob(f'{celltype}/*/{celltype}_*_cis')
            ]
        elif common_or_rare == 'rare':
            existing_assoc_results = [
                str(file)
                for file in to_path(results_path).glob(
                    f'{celltype}/*/{celltype}_*_cis_rare.singleAssoc.txt'
                )
            ]
        # concatenates the dataframes using pandas
        results_all_df_list = []
        results_top_snp_df_list = []
        for pv_df in existing_assoc_results:
            df = pd.read_csv(to_path(pv_df), sep='\t')
            # add gene as column before merging
            df['gene'] = (
                pv_df.split('/')[-1]
                .split('.')[0]
                .replace(f'{celltype}_', '')
                .split('_')[0]
            )
            # add SNP info
            df['is_snp'] = [
                len(df['Allele1'].values[i]) == 1 and len(df['Allele2'].values[i]) == 1
                for i in range(df.shape[0])
            ]
            results_all_df_list.append(df)
            # select top SNP
            df_snp = df[df['is_snp']]
            df_top_snp = df_snp[df_snp['p.value'] == df_snp['p.value'].min()]
            results_top_snp_df_list.append(df_top_snp)
        results_all_df = pd.concat(results_all_df_list)
        results_top_snp_df = pd.concat(results_top_snp_df_list)

        # save
        results_all_file = f'{results_path}/summary_stats/{celltype}_{common_or_rare}_all_cis_raw_pvalues.tsv'
        results_all_df.to_csv(results_all_file, sep='\t')

        results_top_snp_file = f'{results_path}/summary_stats/{celltype}_{common_or_rare}_top_snp_cis_raw_pvalues.tsv'
        results_top_snp_df.to_csv(results_top_snp_file, sep='\t')

        # plot histograms
        # all results
        plt.hist(results_all_df['p.value'])
        plt.savefig('histo.png')
        gcs_path_p = output_path(
            f'plots/pvalues_histo/{celltype}_{common_or_rare}_{title}_all.png',
            'analysis',
        )
        hl.hadoop_copy('histo.png', gcs_path_p)

        # QQ plots
        expected_pvals_all = np.random.uniform(
            low=0, high=1, size=results_all_df.shape[0]
        )
        # all results
        x_all = -np.log10(np.sort(expected_pvals_all))
        y_all = -np.log10(np.sort(results_all_df['p.value']))

        plt.figure(figsize=(8, 8))
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(x_all, y_all)
        fig.savefig('qqplot.png')
        gcs_path_p = output_path(
            f'plots/pvalues_qqplot/{celltype}_{common_or_rare}_{title}_all.png',
            'analysis',
        )
        hl.hadoop_copy('qqplot.png', gcs_path_p)


if __name__ == '__main__':
    plot_pvalues()
