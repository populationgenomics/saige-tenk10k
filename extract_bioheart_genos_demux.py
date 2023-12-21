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

DEFAULT_JOINT_CALL_MT = dataset_path(
    'vds/0c401a0c6aee07549db7d85f88698fbaae4183_990.vds'
)
HAIL_IMAGE = get_config()['images']['scanpy']

# read hail matrix table object (WGS data)
init_batch()

mt_path = DEFAULT_JOINT_CALL_MT
cv_demux_vcf_path = output_path('demux_vcf_common_variants')

mt = hl.read_matrix_table(mt_path)
logging.info(f'Number of total loci: {mt.count()[0]}')
logging.info(f'Number of total samples: {mt.count()[1]}')

mt = mt.filter_rows(
    (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
    & (hl.len(mt.alleles) == 2)  # remove hom-ref
    & (mt.n_unsplit_alleles == 2)  # biallelic (exclude multiallelic)
    & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
)

mt = hl.variant_qc(mt)

# common variants only
cv_mt = mt.filter_rows(mt.variant_qc.AF[0] > 0.05)

# export to plink common variants only for sparse GRM
export_vcf(cv_mt, cv_demux_vcf_path, ind_id=cv_mt.s)
