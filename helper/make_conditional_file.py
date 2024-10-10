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
    python3 helper/make_conditional_file.py \
        --celltypes='B_naive' \
        --results-path='gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files' \
        --conditional-files-output-path='gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files/conditioning_on_top_one_variant_from_cv_test_round1'
"""

import click

import numpy as np
import pandas as pd
from cpg_utils.hail_batch import init_batch

from multipy.fdr import qvalue

# if image doesn't work, copy function from github
# https://github.com/puolival/multipy/blob/master/multipy/fdr.py
from scipy.interpolate import UnivariateSpline


# def qvalue(pvals, threshold=0.05, verbose=True):
#     """Function for estimating q-values from p-values using the Storey-
#     Tibshirani q-value method (2003).

#     Input arguments:
#     ================
#     pvals       - P-values corresponding to a family of hypotheses.
#     threshold   - Threshold for deciding which q-values are significant.

#     Output arguments:
#     =================
#     significant - An array of flags indicating which p-values are significant.
#     qvals       - Q-values corresponding to the p-values.
#     """

#     """Count the p-values. Find indices for sorting the p-values into
#     ascending order and for reversing the order back to original."""
#     m, pvals = len(pvals), np.asarray(pvals)
#     ind = np.argsort(pvals)
#     rev_ind = np.argsort(ind)
#     pvals = pvals[ind]

#     # Estimate proportion of features that are truly null.
#     kappa = np.arange(0, 0.96, 0.01)
#     pik = [sum(pvals > k) / (m * (1 - k)) for k in kappa]
#     cs = UnivariateSpline(kappa, pik, k=3, s=None, ext=0)
#     pi0 = float(cs(1.0))
#     if verbose:
#         print('The estimated proportion of truly null features is %.3f' % pi0)

#     """The smoothing step can sometimes converge outside the interval [0, 1].
#     This was noted in the published literature at least by Reiss and
#     colleagues [4]. There are at least two approaches one could use to
#     attempt to fix the issue:
#     (1) Set the estimate to 1 if it is outside the interval, which is the
#         assumption in the classic FDR method.
#     (2) Assume that if pi0 > 1, it was overestimated, and if pi0 < 0, it
#         was underestimated. Set to 0 or 1 depending on which case occurs.

#     Here we have chosen the first option, since it is the more conservative
#     one of the two.
#     """
#     if pi0 < 0 or pi0 > 1:
#         pi0 = 1
#         print('Smoothing estimator did not converge in [0, 1]')

#     # Compute the q-values.
#     qvals = np.zeros(np.shape(pvals))
#     qvals[-1] = pi0 * pvals[-1]
#     for i in np.arange(m - 2, -1, -1):
#         qvals[i] = min(pi0 * m * pvals[i] / float(i + 1), qvals[i + 1])

#     # Test which p-values are significant.
#     significant = np.zeros(np.shape(pvals), dtype='bool')
#     significant[ind] = qvals < threshold

#     """Order the q-values according to the original order of the p-values."""
#     qvals = qvals[rev_ind]
#     return significant, qvals


@click.command()
@click.option('--celltypes', required=True, help='separated by comma')
@click.option('--results-path', required=True)
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
            # subset to genes in current significant_gene_df
            conditional_df = old_conditional_df[
                old_conditional_df['gene'].isin(significant_gene_df['gene'])
            ]
            # add top marker from current significant_gene_df to the condition in those
            for gene in old_conditional_df['gene']:
                condition = conditional_df[conditional_df['gene'] == gene][
                    'variants_to_condition_on'
                ].split(',')
                new_top_variant = conditional_df[conditional_df['gene'] == gene][
                    'top_MarkerID'
                ]
                # add and reorder by genomic coordinates
                condition.append(new_top_variant)
                conditional_df.update(
                    pd.Series(
                        [','.join(sorted(condition))], name='variants_to_condition_on'
                    )
                )
        # write conditional_df to file
        conditional_file = (
            f'{conditional_files_output_path}/{celltype}_conditional_file.tsv'
        )
        conditional_df.to_csv(conditional_file, sep='\t')


if __name__ == '__main__':
    make_condition_file()
