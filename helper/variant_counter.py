#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter,import-error,no-name-in-module
# ruff: noqa: PLR2004

"""
This script aims to count the total number of variants
from a vds object provided (e.g. bioheart, tob-wgs or both)

similar to tob-wgs/scripts/rv_expression_association/count_variants.py

To run:

analysis-runner \
    --description "count common, low freq and rare variant VCFs" \
    --dataset "bioheart" \
    --access-level "full" \
    --output-dir "saige-qtl/" \
    python3 variant_counter.py --vds-path=gs://cpg-bioheart-main/ashg2024/tenk10k1-0_qc_pass.vds \
        --donors-to-keep=gs://cpg-tob-wgs-main-analysis/large_cohort/tob-wgs1-0/sample_qc.ht/ \
        --donors-to-exclude=gs://cpg-bioheart-main/large_cohort/tenk10k1-1/relateds_to_drop.ht
"""

import click

import hail as hl
import pandas as pd

from cpg_utils.hail_batch import init_batch, output_path


@click.command()
@click.option('--vds-path', required=True)
@click.option('--donors-to-keep', default='all')
@click.option('--donors-to-exclude', default='none')
@click.option('--cv-maf-threshold', default=0.05)
@click.option('--rv-maf-threshold', default=0.01)
@click.option('--exclude-multiallelic', default=False)
@click.option('--exclude-indels', default=False)
def count_variants(
    vds_path: str,
    donors_to_keep: str,
    donors_to_exclude: str,
    cv_maf_threshold: float,
    rv_maf_threshold: float,
    exclude_multiallelic: bool,
    exclude_indels: bool,
):
    """
    reads the VDS, converts to MT,
    if set to true exclude indels and multiallelic snps
    and prints the number of remaining variants
    split between common, low-frequency and rare at given thresholds
    """
    # read VDS object (WGS data)
    init_batch()
    vds_name = vds_path.split('/')[-1].split('.')[0]
    print(f'Counting variants for {vds_name}')
    vds = hl.vds.read_vds(vds_path)
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    # convert to hail matrix table
    mt = hl.vds.to_dense_mt(vds)

    if donors_to_keep != 'all':
        keep_samples_table = hl.read_table(donors_to_keep)
        keep_samples = hl.literal(keep_samples_table.s.collect())
        mt = mt.filter_cols(keep_samples.contains(mt['s']))

    if donors_to_exclude != 'none':
        exclude_samples_table = hl.read_table(donors_to_exclude)
        exclude_samples = hl.literal(exclude_samples_table.s.collect())
        mt = mt.filter_cols(~exclude_samples.contains(mt['s']))

    # filter out loci & variant QC
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # remove hom-ref
    if exclude_multiallelic:  # biallelic only (exclude multiallelic)
        print(f'Number of variants before excluding multiallelic loci: {mt.count()[0]}')
        mt = mt.filter_rows(~(mt.was_split))
        print(f'Number of variants after excluding multiallelic loci: {mt.count()[0]}')
    if exclude_indels:  # SNPs only (exclude indels)
        print(f'Number of variants before excluding indels: {mt.count()[0]}')
        mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
        print(f'Number of variants after excluding indels: {mt.count()[0]}')

    # compute allele frequencies as part of variant qc
    mt = hl.variant_qc(mt)

    # select common, low-frequency and rare variants
    cv_mt = mt.filter_rows(hl.min(mt.variant_qc.AF) >= cv_maf_threshold)
    lf_mt = mt.filter_rows(
        (hl.min(mt.variant_qc.AF) >= rv_maf_threshold)
        & (hl.min(mt.variant_qc.AF) < cv_maf_threshold)
    )
    rv_mt = mt.filter_rows(hl.min(mt.variant_qc.AF) < rv_maf_threshold)

    # count up both donors and variants
    n_common_vars, n_donors = cv_mt.count()
    n_low_frequency_vars = lf_mt.count()[0]
    n_rare_vars = rv_mt.count()[0]

    print(f'Donor count: {n_donors}')

    print(f'Common variant (MAF>={cv_maf_threshold}) count: {n_common_vars}')
    print(
        f'low-frequency variant (MAF >={rv_maf_threshold} and <{cv_maf_threshold}) count: {n_low_frequency_vars}'
    )
    print(f'Rare variant (MAF<{rv_maf_threshold}) count: {n_rare_vars}')

    variant_counter_df = pd.DataFrame(
        [
            {
                'vds_name': vds_name,
                'donor_count': n_donors,
                'rare_variant_maf_threshold': rv_maf_threshold,
                'common_variant_maf_threshold': cv_maf_threshold,
                f'rare_variant_count (MAF<{rv_maf_threshold})': n_rare_vars,
                f'low_frequency_variant_count (MAF >={rv_maf_threshold} and <{cv_maf_threshold})': n_low_frequency_vars,
                f'common_variant_count (MAF>={cv_maf_threshold})': n_common_vars,
            }
        ]
    )
    # save variant counts to file
    variant = 'variant'
    if exclude_multiallelic:
        variant = f'no_multiallelic_{variant}'
    if exclude_indels:
        variant = f'no_indels_{variant}'
    variant_counter_out_file = output_path(
        f'{vds_name}_mafs_{rv_maf_threshold}_{cv_maf_threshold}_{variant}_counts.csv',
        'analysis',
    )
    variant_counter_df.to_csv(variant_counter_out_file)


if __name__ == '__main__':
    count_variants()
