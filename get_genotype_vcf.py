#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script will

- extract common variants
- export as VCF (one per chrom)

this will be used as input for the
SAIGE-QTL association pipeline.
"""

# python modules
from cpg_utils.hail_batch import (
    dataset_path,
    get_config,
    init_batch,
    output_path,
)

import hail as hl
from hail.methods import export_vcf


BIOHEART_TOB_VDS = dataset_path('.vds')
HAIL_IMAGE = get_config()['images']['scanpy']
CHROMOSOMES = 'chr22'

# read hail matrix table object (WGS data)
init_batch()

vds_path = BIOHEART_TOB_VDS

vds = hl.vds.read_vds(vds_path)

for chromosome in CHROMOSOMES:

    # consider only relevant chromosome
    vds = hl.vds.filter_chromosomes(vds, keep=chromosome)

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
    cv_mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.05)

    # remove fields not in the VCF
    cv_mt = cv_mt.drop('gvcf_info')

    # export to plink common variants only
    cv_vcf_path = output_path(f'{chromosome}_common_variants.vcf.bgz')
    export_vcf(cv_mt, cv_vcf_path)

    # figure out how to create index file (.csi)
