#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script will

- extract common, biallelic SNPs
- export as VCF

this excludes non-variable loci since they will not help
in differentiating individuals, multi-allelic SNPs, as they
are more likely to be artefact and harder to find in RNA data,
and indels which are also not useful for demultiplexing

this will be used to demultiplex scRNA-seq
data so that we can map cells to donors.

To run:

analysis-runner \
    --description "get common variant vcf for demultiplexing" \
    --dataset "bioheart" \
    --access-level "standard" \
    --output-dir "saige-qtl/demux_files/" \
    python3 helper/extract_bioheart_genos_demux.py

"""

import click

from cpg_utils.hail_batch import (
    dataset_path,
    init_batch,
    output_path,
)

import hail as hl
from hail.methods import export_vcf

import logging
logging.basicConfig(level=logging.INFO)

@click.command()
@click.option('--vds-path', default=dataset_path('vds/bioheart1-0.vds'))
@click.option(
    '--cv-demux-vcf-path',
    default=output_path('bioheart_demux_vcf_common_variants.vcf.bgz'),
)
@click.option('--cv-maf-threshold', default=0.05)
def get_demux_vcf(
    vds_path: str,
    cv_demux_vcf_path: str,
    cv_maf_threshold: float,
):
    """
    reads the VDS, extract common SNPs,
    writes to VCF
    """
    init_batch()

    vds = hl.vds.read_vds(vds_path)

    # split multiallelic loci
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    # densify to matrix table object
    mt = hl.vds.to_dense_mt(vds)

    # filter out loci & variant QC
    mt = mt.filter_rows(
        ~(mt.was_split)  # biallelic (exclude multiallelic)
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
    )
    mt = hl.variant_qc(mt)

    # common variants only
    cv_mt = mt.filter_rows(hl.min(mt.variant_qc.AF) > cv_maf_threshold)

    # remove fields not in the VCF
    cv_mt = cv_mt.drop('gvcf_info')

    # export to plink common variants only for demultiplexing
    export_vcf(cv_mt, cv_demux_vcf_path)


if __name__ == '__main__':
    get_demux_vcf()
