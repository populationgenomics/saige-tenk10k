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

__author__ = 'annacuomo'

BIOHEART_JOINT_CALL_VDS = dataset_path(
    'vds/0c401a0c6aee07549db7d85f88698fbaae4183_990.vds'
)
HAIL_IMAGE = get_config()['images']['scanpy']

# read hail matrix table object (WGS data)
init_batch()

vds_path = BIOHEART_JOINT_CALL_VDS

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
cv_mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.05)

# remove fields not in the VCF
cv_mt = cv_mt.drop('gvcf_info')

# export to plink common variants only for demultiplexing
cv_demux_vcf_path = output_path('demux_vcf_common_variants.vcf.bgz')
export_vcf(cv_mt, cv_demux_vcf_path)
