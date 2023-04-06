#!/usr/bin/env python3
# pylint: disable=broad-exception-raised,import-error,import-outside-toplevel,missing-module-docstring,no-value-for-parameter,too-many-arguments,too-many-branches,too-many-locals,too-many-statements,wrong-import-order,wrong-import-position

__author__ = 'annacuomo'

"""
Hail Batch workflow for the rare-variant association analysis, including:

- perform sample and variant QC
- get relevant variants around a gene and export genotypes as plink files,
- generate other input files for association tests (phenotype, covariates, groups),
- run association tests.

More details in README
output files in tob_wgs_genetics/saige_qtl/output
"""

# import python modules
import sys

import logging

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    copy_common_env,
    dataset_path,
    get_config,
    init_batch,
    output_path,
    remote_tmpdir,
)

import click
import pandas as pd

import hail as hl
import hailtop.batch as hb

from sample_metadata.apis import SequenceApi

# use logging to print statements, display at info level
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr,
)

seqapi = SequenceApi()

DEFAULT_JOINT_CALL_MT = dataset_path('mt/v7.mt')

HAIL_IMAGE = get_config()['workflow']['driver_image']


# region SUBSET_SAMPLES


def remove_sc_outliers(df, outliers=None):
    """
    Remove outlier samples, as identified by single-cell analysis
    """
    if outliers is None:
        outliers = ['966_967', '88_88']
    else:
        outliers = outliers.extend(['966_967', '88_88'])
    df = df[-df['OneK1K_ID'].isin(outliers)]

    return df


# extract bone marrow samples
# it's a sequence metadata vs sample??
def get_bone_marrow_samples():
    """
    Extract TOB bone marrow samples (vs PBMCs)
    """
    bm_samples = seqapi.get_samples(
        body_get_samples={
            'seq_meta': {'Primary study': 'Pilot/bone marrow'},
            'projects': ['tob-wgs'],
        }
    )
    return bm_samples


# remove duplicated samples based on TOB IDs
# CPG4994, CPG67264 both the same individual (TOB1282)
# CPG5066, CPG67504 both TOB1289
# in both cases keep the latter which is the resequenced version
# looking for better than manual extraction of these two
def get_duplicated_samples():
    """
    Extract duplicated samples for same individual
    """
    duplicated_samples = ['CPG4994', 'CPG5066']
    return duplicated_samples


def get_non_tob_samples():
    """
    Extract outsider samples not from this cohort
    """
    outsiders = ['NA12878', 'NA12891', 'NA12892', 'Syndip']
    return outsiders


# endregion SUBSET_SAMPLES

# region SUBSET_VARIANTS


# only needs to be run once for a given cohort (e.g., OneK1K)
def filter_variants(
    mt_path: str,  # 'mt/v7.mt'
    samples: list[str],
    output_rv_mt_path: str,  # 'tob_wgs/densified_rv_only.mt'
    output_cv_mt_path: str,  # 'tob_wgs/densified_cv_only.mt'
    vre_plink_path: str,  # 'tob_wgs/vr_plink_2000_variants
    cv_maf_threshold: float = 0.01,
    rv_maf_threshold: float = 0.05,
    vre_mac_threshold: int = 20,
    vre_n_markers: int = 2000,
):
    """Subset hail matrix table

    Input:
    - joint call hail matrix table
    - set of samples for which we have scRNA-seq data
    - file paths for outputs
    - MAF thresholds to define common / rare variants

    Output 1&2:
    subset hail matrix table, containing only variants that:
    i) are not ref-only, ii) biallelic, iii) meet QC filters,
    and samples that are contained in sc sample list.

    Then, in output1 variants are also rare (freq < rv_maf_threshold)
    and in output2 they are common (freq > cv_maf_threshold)

    Output 3:
    plink file containing a random subset of 2,000 variants that satisfy i),ii),iii)
    that are additionally sufficiently common (MAC>20) and not in LD
    """
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # subset to relevant samples (samples we have scRNA-seq data for)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # densify
    mt = hl.experimental.densify(mt)

    # add column filters
    bm_samples = get_bone_marrow_samples()
    dup_samples = get_duplicated_samples()
    out_samples = get_non_tob_samples()
    # qc_samples = get_qced_out_samples()
    mt = mt.filter_cols(mt.s in bm_samples)
    mt = mt.filter_cols(mt.s in dup_samples)  # merge with above?
    mt = mt.filter_cols(mt.s in out_samples)

    # filter out low quality variants and consider biallelic SNPs only
    # (no multi-allelic, no ref-only, no indels)
    mt = mt.filter_rows(  # check these filters!
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic (revisit)
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (revisit)
    )

    mt = hl.variant_qc(mt)

    # subset variants for variance ratio estimation
    # minor allele count (MAC) > 20
    tot_counts = mt.variant_qc.AC.sum()
    vre_mt = mt.filter_rows(
        (mt.variant_qc.AC[1] > vre_mac_threshold) & (mt.variant_qc.AC[1] < tot_counts)
        | (mt.variant_qc.AF[1] < (tot_counts - vre_mac_threshold))
        & (mt.variant_qc.AC[1] > 0)
    )
    # perform LD pruning
    vre_mt = vre_mt.sample_rows(
        p=0.01
    )  # in case this is very costly, subset first a bit
    pruned_variant_table = hl.ld_prune(vre_mt.GT, r2=0.2, bp_window_size=500000)
    vre_mt = vre_mt.filter_rows(hl.is_defined(pruned_variant_table[vre_mt.row_key]))
    # randomly sample {vre_n_markers} variants
    vre_mt = vre_mt.sample_rows((vre_n_markers * 1.1) / vre_mt.count[0])
    vre_mt = vre_mt.head(vre_n_markers)

    # export to plink common variants only for sparse GRM
    from hail.methods import export_plink

    export_plink(vre_mt, vre_plink_path, ind_id=vre_mt.s)

    # filter common variants for single-variant association
    cv_mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > cv_maf_threshold) & (mt.variant_qc.AF[1] < 1)
        | (mt.variant_qc.AF[1] < (1 - cv_maf_threshold)) & (mt.variant_qc.AF[1] > 0)
    )
    cv_mt.write(output_cv_mt_path, overwrite=True)
    logging.info(
        f'No common (freq>{cv_maf_threshold}), biallelic SNPs: {cv_mt.count()[0]}'
    )

    # filter rare variants only (MAF < 5%)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] < rv_maf_threshold) & (mt.variant_qc.AF[1] > 0)
        | (mt.variant_qc.AF[1] > (1 - rv_maf_threshold)) & (mt.variant_qc.AF[1] < 1)
    )
    mt.write(output_rv_mt_path, overwrite=True)
    logging.info(f'No rare (freq<{rv_maf_threshold}), biallelic SNPs: {mt.count()[0]}')


# endregion SUBSET_VARIANTS


config = get_config()


@click.command()
@click.option(
    '--sample-mapping-file-tsv',
    default='scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv',  # to be updated
)
@click.option('--mt-path', default=DEFAULT_JOINT_CALL_MT)
def saige_pipeline(
    sample_mapping_file_tsv: str,
    mt_path: str,
):
    """
    Run entire pipeline
    """
    sb = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch('SAIGE-QTL pipeline', backend=sb)

    # extract individuals for which we have single-cell (sc) data
    sample_mapping_file = pd.read_csv(dataset_path(sample_mapping_file_tsv), sep='\t')
    # we may want to exclude these from the smf directly
    sample_mapping_file = remove_sc_outliers(sample_mapping_file)
    # check column names - CPG_ID would be better?
    sc_samples = sample_mapping_file['InternalID'].unique()

    # filter to QC-passing, rare, biallelic variants
    output_rv_mt_path = output_path('densified_rv_only.mt')
    output_cv_mt_path = output_path('densified_cv_only.mt')
    vre_plink_path = output_path('vr_plink_20k_variants')
    if not to_path(output_rv_mt_path).exists():

        filter_job = batch.new_python_job(name='MT filter job')
        copy_common_env(filter_job)
        filter_job.image(HAIL_IMAGE)
        filter_job.call(
            filter_variants,
            mt_path=mt_path,
            samples=list(sc_samples),
            output_rv_mt_path=output_rv_mt_path,
            output_cv_mt_path=output_cv_mt_path,
            vre_plink_path=vre_plink_path,
        )

    else:
        logging.info('File already exists no need to filter')
        filter_job = None

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    saige_pipeline()
