#!/usr/bin/env python3
"""
This script summarise SAIGE-QTL raw p-value results
and plots QQ plots and histograms of the p-values

To run:

analysis-runner \
    --description "plot p-values" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/" \
    python3 summarise_and_qq_plotter.py \
        --celltype='B_naive' \
        --results-path=gs://cpg-bioheart-main-analysis/saige-qtl/bioheart_n990_and_tob_n1055/output_files/sample_perm0/output_files \
        --title='Shuffled p-values (SAIGE-QTL pipeline)'
"""

import click
import matplotlib.pyplot as plt
import hail as hl

# import numpy as np
import pandas as pd
from cpg_utils import to_path
from cpg_utils.hail_batch import output_path

# from cpg_utils.hail_batch import init_batch, output_path
# from bokeh.plotting import output_file, save


@click.command()
@click.option('--celltype', required=True)
@click.option('--results-path', required=True)
# @click.option('--title', default='SAIGE-QTL pipeline p-values')
def plot_pvalues(
    celltype: str,
    results_path: str,
    # title: str,
):
    """
    combines the results for a given cell type,
    saving both all results and top SNP per gene results
    plots both a histogram and a QQ plot of the association p-values
    """
    # init_batch()

    # collect all raw p-value files
    existing_cv_assoc_results = [
        str(file) for file in to_path(results_path).glob(f'{celltype}_*_cis')
    ]
    # concatenates the dataframes using pandas
    results_all_df_list = []
    results_top_snp_df_list = []
    for pv_df in existing_cv_assoc_results:
        df = pd.read_csv(to_path(pv_df), sep='\t')
        # add gene as column before merging
        df['gene'] = (
            pv_df.replace(f'{results_path}/', '')
            .replace(f'{celltype}_', '')
            .replace('_cis', '')
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
    results_all_file = (
        f'{results_path}/summary_stats/{celltype}_all_cis_raw_pvalues.tsv'
    )
    results_all_df.to_csv(results_all_file, sep='\t')

    results_top_snp_file = (
        f'{results_path}/summary_stats/{celltype}_top_snp_cis_raw_pvalues.tsv'
    )
    results_top_snp_df.to_csv(results_top_snp_file, sep='\t')

    # plot histograms
    p_hist_all = plt.hist(results_all_df['p.value'])
    # p_hist_top = plt.hist(results_top_snp_df['p.value'])

    # # QQ plots
    # expected_pvals_all = np.random.uniform(low=0, high=1, size=results_all_df.shape[0])
    # expected_pvals_top = np.random.uniform(
    #     low=0, high=1, size=results_top_snp_df.shape[0]
    # )

    fig = plt.subplots(figsize=(10, 8))
    fig.save(p_hist_all)
    gcs_path_p = output_path(
        f'plots/pvalues_histo/{celltype}_shuffled.html', 'analysis'
    )
    hl.hadoop_copy('local_histo.html', gcs_path_p)


if __name__ == '__main__':
    plot_pvalues()
