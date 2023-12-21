#!/usr/bin/env python3

__author__ = 'annacuomo'

"""
This script will

- extract common variants
- export as VCF

this will be used to demultiplex scRNA-seq
data so that we can map cells to donors.
"""

import logging

from cpg_utils.hail_batch import (
    dataset_path,
    get_config,
    init_batch,
    output_path,
)

import hail as hl
from hail.methods import export_vcf

BIOHEART_JOINT_CALL_VDS = dataset_path(
    'vds/0c401a0c6aee07549db7d85f88698fbaae4183_990.vds'
)
HAIL_IMAGE = get_config()['images']['scanpy']

# read hail matrix table object (WGS data)
init_batch()

vds_path = BIOHEART_JOINT_CALL_VDS
cv_demux_vcf_path = output_path('demux_vcf_common_variants')

vds = hl.read_vds(vds_path)
logging.info(f'Number of total loci: {vds.count()[0]}')
logging.info(f'Number of total samples: {vds.count()[1]}')

vds = vds.filter_rows(
    (hl.len(hl.or_else(vds.filters, hl.empty_set(hl.tstr))) == 0)  # QC
    & (hl.len(vds.alleles) == 2)  # remove hom-ref
    & (vds.n_unsplit_alleles == 2)  # biallelic (exclude multiallelic)
    & (hl.is_snp(vds.alleles[0], vds.alleles[1]))  # SNPs (exclude indels)
)

vds = hl.variant_qc(vds)

# common variants only
cv_vds = vds.filter_rows(vds.variant_qc.AF[0] > 0.05)

# export to plink common variants only for sparse GRM
export_vcf(cv_vds, cv_demux_vcf_path, ind_id=cv_vds.s)
