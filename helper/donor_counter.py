#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter,import-error,no-name-in-module
# ruff: noqa: PLR2004

"""
This script aims to count the number of donors present in the different
input files for SAIGE-QTL to identify why they do not match and hence
the run is failing (pt 2): https://batch.hail.populationgenomics.org.au/batches/478962

To run:

analysis-runner \
    --description "count VCF / VRE donors" \
    --dataset "bioheart" \
    --access-level "full" \
    --output-dir "saige-qtl/" \
    python3 donor_counter.py --chrom=chr22 \
    --genos-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/genotypes/v2/vds-tenk10k1-0
"""

import click

import hail as hl

# from cpg_utils.hail_batch import init_batch


@click.command()
@click.option('--chrom', required=True)
@click.option('--genos-path', required=True)
def count_variants(
    chrom: str,
    genos_path: str,
):
    """
    reads the VCF files (common & rare) as MTs,
    reads the VRE .fam file,
    extracts list of donors from all files, compares
    """
    # init_batch()
    vcf_path_common = f'{genos_path}/{chrom}_common_variants.vcf.bgz'
    mt_common = hl.import_vcf(vcf_path_common)

    # vcf_path_rare = f'{genos_path}/{chrom}_rare_variants.vcf.bgz'
    # mt_rare = hl.import_vcf(vcf_path_rare)

    mt_plink = hl.import_plink(
        bed=f'{genos_path}/vre_plink_2000_variants.bed',
        bim=f'{genos_path}/vre_plink_2000_variants.bim',
        fam=f'{genos_path}/vre_plink_2000_variants.fam',
        reference_genome='GRCh38',
    )

    # extract donors
    mt_common_donors = mt_common.s.collect()
    # mt_rare_donors = mt_rare.s.collect()
    mt_plink_donors = mt_plink.s.collect()

    # count donors
    print(f'Donor count (common variant vcf): {len(mt_common_donors)}')
    # print(f'Donor count (rare variant vcf): {len(mt_rare_donors)}')
    print(f'Donor count (VRE plink file): {len(mt_plink_donors)}')

    # compare lists
    # print(
    #     f'Do {vcf_path_common} and {vcf_path_rare} contain the same donors? {mt_common_donors.sort() == mt_rare_donors.sort()}'
    # )
    print(
        f'Do {vcf_path_common} and {genos_path}/vre_plink_2000_variants.fam contain the same donors? {mt_common_donors.sort() == mt_plink_donors.sort()}'
    )
    # print(
    #     f'Do {vcf_path_rare} and {genos_path}/vre_plink_2000_variants.fam contain the same donors? {mt_rare_donors.sort() == mt_plink_donors.sort()}'
    # )


if __name__ == '__main__':
    count_variants()
