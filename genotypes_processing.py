#!/usr/bin/env python3
# pylint: disable=broad-exception-raised,import-error,import-outside-toplevel,missing-module-docstring,no-value-for-parameter,too-many-arguments,too-many-branches,too-many-locals,too-many-statements,wrong-import-order,wrong-import-position

__author__ = 'annacuomo'

"""
Hail Batch workflow to extract relevant variants to test.
This script will:

- perform sample QC
- perform variant QC
- export filtered matrix table
- export genotype subset as plink files,

More details in README
output files in tob_wgs_genetics/saige_qtl/input
"""

# import python modules
import sys

import logging

import random

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    get_config,
    init_batch,
    output_path,
)

import click
import pandas as pd

import hail as hl


from metamist.apis import ParticipantApi, SequencingGroupApi


# use logging to print statements, display at info level
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr,
)

papi = ParticipantApi()


DEFAULT_JOINT_CALL_MT = dataset_path('mt/v7.mt')
HAIL_IMAGE = get_config()['workflow']['driver_image']


# region SAMPLES_SUBSETTING_FUNCTIONS


def remove_sc_outliers(df, outliers=None):
    """
    Remove outlier samples, as identified by single-cell analysis
    """
    if outliers is None:
        # by default, for TOB, we remove two samples
        # 966_967 for extremely low cell number
        # 88_88 for abnormal B cell composition
        outliers = ['966_967', '88_88']
    else:
        outliers = outliers.extend(['966_967', '88_88'])
    df = df[-df['OneK1K_ID'].isin(outliers)]

    return df


# for the 3 functions below, not very confident on the use of
# metamist, so I have retained the manual selection of
# specific samples (commented out)
# as a note: all of this should be done within the processing
# pipeline at some stage, so perhaps does not need to be perfect


def get_bone_marrow_sequencing_groups():
    """
    Extract TOB bone marrow samples (vs PBMCs)
    """
    from metamist.graphql import gql, query

    _query = gql(
        """
    query MyQuery {
        project(name: "tob-wgs") {
            sequencingGroups {
                id
                meta
                assays {
                    id
                    meta
                }
            }
        }
    }
    """
    )
    logging.getLogger().setLevel(logging.WARN)
    sequencing_groups = query(_query)['project']['sequencingGroups']
    logging.getLogger().setLevel(logging.INFO)

    bm_sequencing_groups = []
    for sg in sequencing_groups:
        sg_id = sg['id']
        for assay in sg['assays']:
            if assay['meta'].get('Primary study') == 'Pilot/bone marrow':
                bm_sequencing_groups.append(sg_id)
                continue
    logging.info(f'Number of bone marrow samples: {len(set(bm_sequencing_groups))}')
    return set(bm_sequencing_groups)


# remove duplicated samples based on TOB IDs
# CPG4994, CPG67264 both the same individual (TOB1282)
# CPG5066, CPG67504 both TOB1289
# in both cases keep the latter which is the resequenced version
def get_duplicated_samples(mt: hl.MatrixTable) -> set:
    """
    Extract duplicated samples for same individual
    i.e., keep "-PBMC" version (now keeping most recent)
    """
    sgapi = SequencingGroupApi()
    sample_sg_map = sgapi.get_all_sequencing_group_ids_by_sample_by_type(
        project='tob-wgs'
    )
    sams = papi.get_external_participant_id_to_internal_sample_id(project='tob-wgs')
    # Keep the most recent sample for each participant in the project
    latest_samples = set(dict(sams).values())

    # Get the SG ID associated with each of the most recent samples
    keep = [sample_sg_map[sample_id] for sample_id in latest_samples]
    keep = [list(sg.values())[0] for sg in keep]
    # lint complains about syntax below even though it works?
    keep = set([sublist for list in keep for sublist in list])

    matrix_samples = set(mt.s.collect())
    dup_samples = matrix_samples.difference(keep)
    logging.info(f'Number of duplicated samples: {len(set(dup_samples))}')
    print(set(dup_samples))
    # if set(dup_samples) != {'CPG4994', 'CPG5066'}:
    #     logging.info("Not the right samples, check this function")
    #     return set()
    # return {'CPG4994', 'CPG5066'}
    return set(dup_samples)


def get_non_tob_samples(mt: hl.MatrixTable) -> set:
    """
    Extract outsider samples not from this cohort
    (only included for comparison purpose)
    """
    sgapi = SequencingGroupApi()
    # Return Sample IDs mapped to seq type, sequencing group ID
    # (e.g. {'XPG123': {'genome' : ['CPG123']}, {'XPG456': {'exome':['CPG456']}})
    sample_sg_map = sgapi.get_all_sequencing_group_ids_by_sample_by_type(
        project='tob-wgs'
    )
    sgs = [list(sg.values())[0] for sg in sample_sg_map.values()]
    # double-layered list comprehension to flatten
    tob_samples = [sublist for list in sgs for sublist in list]
    matrix_samples = set(mt.s.collect())
    common_samples = set(tob_samples).intersection(matrix_samples)
    logging.info(f'Matrix samples: {matrix_samples}')
    logging.info(f'Common samples: {common_samples}')
    if common_samples == matrix_samples:
        logging.info('No samples to remove, exit')
        return set()
    # non_tob_samples = matrix_samples not in common_samples
    non_tob_samples = matrix_samples.difference(common_samples)
    logging.info(f'Number of non-TOB samples: {len(non_tob_samples)}')
    print(non_tob_samples)
    if non_tob_samples != {'NA12878', 'NA12891', 'NA12892', 'Syndip'}:
        logging.info('Not the right samples, check this function')
        return set()
    # return {'NA12878', 'NA12891', 'NA12892', 'Syndip'}
    return {non_tob_samples}


def get_low_qc_samples(
    metadata_tsv_path='gs://cpg-tob-wgs-test-analysis/joint-calling/v7/meta.tsv',
    contam_rate=0.05,  # defined based on distribution of parameters here
    chimera_rate=0.05,  # as above
):
    """
    Extract samples that did not pass QC
    - high contamination rate
    - high chimera rate
    - ambiguous sex or sex aneuploidy
    - WGS-based QC
    """
    meta = pd.read_csv(metadata_tsv_path, sep='\t')
    samples_contam = set(meta[meta['r_contamination'] > contam_rate]['s'])
    logging.info(
        f'Number of samples w. contam. rate > {contam_rate}: {len(samples_contam)}'
    )
    samples_chim = set(meta[meta['r_chimera'] > chimera_rate]['s'])
    logging.info(
        f'Number of samples with chimera rate > {chimera_rate}: {len(samples_chim)}'
    )
    samples_sex = {
        *set(meta[meta['hard_filters'] == 'ambiguous_sex']['s']),
        *set(meta[meta['hard_filters'] == 'sex_aneuploidy']['s']),
    }
    logging.info(
        f'Number of samples with sex aneuploidy or ambiguous sex: {len(samples_sex)}'
    )
    samples_qc = set(meta[meta['qc_metrics_filters'].notnull()]['s'])
    logging.info(f'Number of samples not passing WGS QC: {len(samples_qc)}')
    return {*samples_contam, *samples_chim, *samples_sex, *samples_qc}


# endregion SAMPLES_SUBSETTING_FUNCTIONS

# region VARIANTS_SUBSETTING_FUNCTIONS


# only needs to be run once for a given cohort (e.g., OneK1K / TOB)
def filter_variants(
    mt_path: str,  # 'mt/v7.mt'
    samples_str: str,
    output_mt_path: str,  # 'tob_wgs/filtered.mt'
    vre_plink_path: str,  # 'tob_wgs/vr_plink_2000_variants
    vre_mac_threshold: int = 20,
    vre_n_markers: int = 2000,
):
    """Subset hail matrix table

    Input:
    - joint call hail matrix table
    - set of samples for which we have scRNA-seq data
    - file paths for outputs

    Output 1:
    subset hail matrix table, containing only variants that:
    i) are not ref-only, ii) biallelic, iii) meet QC filters,
    and samples that are contained in sc sample list + QC passing.

    Output 2:
    plink file containing a random subset of 2,000 variants that satisfy i),ii),iii)
    and also post sample QC
    that are additionally sufficiently common (MAC>20) and not in LD
    """
    samples = samples_str.split(',')
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)
    logging.info(f'Number of total loci: {mt.count()[0]}')
    logging.info(f'Number of total samples: {mt.count()[1]}')

    # subset to relevant samples (samples we have scRNA-seq data for)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # densify
    mt = hl.experimental.densify(mt)

    # add sample filters
    bm_samples = get_bone_marrow_sequencing_groups()
    dup_samples = get_duplicated_samples(mt=mt)
    out_samples = get_non_tob_samples(mt=mt)
    qc_samples = get_low_qc_samples()
    filter_samples = {*bm_samples, *dup_samples, *out_samples, *qc_samples}
    # will this work with a set or should it be a list?
    mt = mt.filter_cols(mt.s in filter_samples)
    logging.info(f'Number of samples after filtering: {mt.count()[1]}')
    if mt.count()[1] == 0:
        logging.info('No samples left, exit')
        return

    # filter out low quality variants and consider biallelic SNPs only
    # (no multi-allelic, no ref-only, no indels)
    mt = mt.filter_rows(  # check these filters!
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic (exclude multiallelic)
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
    )

    mt = hl.variant_qc(mt)

    mt.write(output_mt_path, overwrite=True)
    logging.info(f'No QC-passing, biallelic SNPs: {mt.count()[0]}')

    # subset variants for variance ratio estimation
    # minor allele count (MAC) > 20
    vre_mt = mt.filter_rows(mt.variant_qc.AC[0] > vre_mac_threshold)
    # when testing this in a notebook I get 8809223 here I keep getting 0
    logging.info(f'Number of variants post AC filter: {vre_mt.count()[0]}')
    if vre_mt.count()[0] == 0:
        logging.info('No variants left, exit')
        return
    # perform LD pruning
    vre_mt = vre_mt.sample_rows(
        p=0.01
    )  # in case this is very costly, subset first a bit
    logging.info(f'Initial subset of variants: {vre_mt.count()[0]}')
    pruned_variant_table = hl.ld_prune(vre_mt.GT, r2=0.2, bp_window_size=500000)
    vre_mt = vre_mt.filter_rows(hl.is_defined(pruned_variant_table[vre_mt.row_key]))
    logging.info(f'Subset of variants after pruning: {vre_mt.count()[0]}')
    # randomly sample {vre_n_markers} variants
    random.seed(0)
    vre_mt = vre_mt.sample_rows((vre_n_markers * 1.1) / vre_mt.count[0])
    vre_mt = vre_mt.head(vre_n_markers)
    logging.info(f'Subset to {vre_n_markers} variants: {vre_mt.count()[0]}')

    # export to plink common variants only for sparse GRM
    from hail.methods import export_plink

    export_plink(vre_mt, vre_plink_path, ind_id=vre_mt.s)


# endregion VARIANTS_SUBSETTING_FUNCTIONS


@click.command()
@click.option(
    '--sample-mapping-file-tsv',
    default='scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv',  # to be updated
)
@click.option('--mt-path', default=DEFAULT_JOINT_CALL_MT)
def genotypes_pipeline(
    sample_mapping_file_tsv: str,
    mt_path: str,
):
    """
    Run one-off QC filtering pipeline
    """

    # extract individuals for which we have single-cell (sc) data
    sample_mapping_file = pd.read_csv(dataset_path(sample_mapping_file_tsv), sep='\t')
    # we may want to exclude these from the smf directly
    sample_mapping_file = remove_sc_outliers(sample_mapping_file)
    # check column names - CPG_ID would be better?
    sc_samples = ','.join(sample_mapping_file['InternalID'].unique())

    # filter to QC-passing, biallelic SNPs
    output_mt_path = output_path('qc_filtered.mt')
    vre_plink_path = output_path('vr_plink_20k_variants')
    # only exit if both files exist
    if all(to_path(output).exists() for output in [output_mt_path, vre_plink_path]):
        logging.info('File already exists no need to filter')
        return

    filter_variants(
        mt_path=mt_path,
        samples_str=sc_samples,
        output_mt_path=output_mt_path,
        vre_plink_path=vre_plink_path,
    )


if __name__ == '__main__':
    genotypes_pipeline()
