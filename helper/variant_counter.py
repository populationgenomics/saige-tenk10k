#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter,import-error,no-name-in-module
# ruff: noqa: PLR2004

"""
This script aims to count the total number of variants
from a vds object provided (e.g. bioheart, tob-wgs or both)

similar to tob-wgs/scripts/rv_expression_association/count_variants.py

To run:

analysis-runner \
    --description "count common and rare variant VCFs" \
    --dataset "bioheart" \
    --access-level "full" \
    --output-dir "saige-qtl/" \
    python3 variant_counter.py --vds-path=gs://cpg-tob-wgs-main/vds/tob-wgs1-0.vds
"""

import click

import hail as hl

from cpg_utils.hail_batch import init_batch


@click.command()
@click.option('--vds-path', required=True)
@click.option('--cv-maf-threshold', default=0.05)
@click.option('--rv-maf-threshold', default=0.05)
@click.option('--exclude-multiallelic', default=False)
@click.option('--exclude-indels', default=False)
def count_variants(
    vds_path: str,
    cv_maf_threshold: float,
    rv_maf_threshold: float,
    exclude_multiallelic: bool,
    exclude_indels: bool,
):
    """
    reads the VDS, converts to MT,
    if set to true exclude indels and multiallelic snps
    and prints the number of remaining variants
    split between common and rare at given thresholds
    """
    # read VDS object (WGS data)
    init_batch()
    vds_name = vds_path.split('/')[-1].split('.')[0]
    print(f'Counting variants for {vds_name}')
    vds = hl.vds.read_vds(vds_path)
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    # convert to hail matrix table
    mt = hl.vds.to_dense_mt(vds)

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

    # select common and rare variants
    cv_mt = mt.filter_rows(mt.variant_qc.AF[1] >= cv_maf_threshold)
    rv_mt = mt.filter_rows(mt.variant_qc.AF[1] < rv_maf_threshold)

    print(f'Common variant (MAF>={cv_maf_threshold}) count: {cv_mt.count()[0]}')
    print(f'Rare variant (MAF<{rv_maf_threshold}) count: {rv_mt.count()[0]}')


if __name__ == '__main__':
    count_variants()
