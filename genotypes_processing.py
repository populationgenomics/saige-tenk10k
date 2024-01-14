#!/usr/bin/env python3
# pylint: disable=import-error,import-outside-toplevel,no-value-for-parameter

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

import logging
import random
import sys

import click
import hail as hl
import pandas as pd
from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    get_config,
    init_batch,
    output_path,
)
from metamist.apis import ParticipantApi, SequencingGroupApi

__author__ = 'annacuomo'

# use logging to print statements, display at info level
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr,
)

papi = ParticipantApi()


DEFAULT_JOINT_CALL_MT = dataset_path('mt/v7.mt')
HAIL_IMAGE = get_config()['images']['scanpy']


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


def get_bone_marrow_sequencing_groups(mt: hl.MatrixTable) -> set:
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
    sequencing_groups = query(_query)['project']['sequencingGroups']  # type: ignore
    logging.getLogger().setLevel(logging.INFO)

    bm_sequencing_groups = []
    for sg in sequencing_groups:
        sg_id = sg['id']
        for assay in sg['assays']:
            if assay['meta'].get('Primary study') == 'Pilot/bone marrow':
                bm_sequencing_groups.append(sg_id)
                continue
    logging.info(f'Number of bone marrow samples: {len(set(bm_sequencing_groups))}')
    bm_samples = set(bm_sequencing_groups)
    matrix_samples = set(mt.s.collect())
    return bm_samples.intersection(matrix_samples)


# remove duplicated samples based on TOB IDs
# keep the latter which is the resequenced version
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
    keep = set()
    for sample_id in latest_samples:
        keep.update(sample_sg_map[sample_id]['genome'])

    matrix_samples = set(mt.s.collect())
    dup_samples = matrix_samples.difference(keep)
    logging.info(f'Number of duplicated samples: {len(set(dup_samples))}')
    print(set(dup_samples))
    return set(dup_samples)


def get_non_tob_samples(mt: hl.MatrixTable) -> set:
    """
    Extract outsider samples not from this cohort
    (only included for comparison purpose)
    """
    sgapi = SequencingGroupApi()
    # Return Sample IDs mapped to seq type, sequencing group ID
    sample_sg_map = sgapi.get_all_sequencing_group_ids_by_sample_by_type(
        project='tob-wgs'
    )
    tob_samples = set()
    for sample_id in sample_sg_map:
        tob_samples.update(sample_sg_map[sample_id]['genome'])
    matrix_samples = set(mt.s.collect())
    logging.info(f'Matrix samples: {matrix_samples}')
    common_samples = set(tob_samples).intersection(matrix_samples)
    logging.info(f'Common samples: {common_samples}')
    if common_samples == matrix_samples:
        logging.info('No samples to remove, exit')
        return set()
    non_tob_samples = matrix_samples.difference(common_samples)
    logging.info(f'Number of non-TOB samples: {len(non_tob_samples)}')
    print(non_tob_samples)
    if non_tob_samples != {'NA12878', 'NA12891', 'NA12892', 'syndip'}:
        logging.info('Not the right samples, check non-TOB function')
        return set()
    return non_tob_samples


# cutoffs currently based on:
# https://github.com/populationgenomics/joint-calling/blob/main/joint_calling/filter_cutoffs.yaml
def get_low_qc_samples(
    mt: hl.MatrixTable,
    metadata_tsv='joint-calling/v7/meta.tsv',
    contam_rate=0.04,  # defined based on distribution of parameters in TOB
    chimera_rate=0.05,  # as above
) -> set:
    """
    Extract samples that did not pass QC
    - high contamination rate
    - high chimera rate
    - ambiguous sex or sex aneuploidy
    - WGS-based QC
    """
    metadata_tsv_path = dataset_path(metadata_tsv, 'analysis')
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
    low_qc_samples = {*samples_contam, *samples_chim, *samples_sex, *samples_qc}
    matrix_samples = set(mt.s.collect())
    return low_qc_samples.intersection(matrix_samples)


# endregion SAMPLES_SUBSETTING_FUNCTIONS

# region VARIANTS_SUBSETTING_FUNCTIONS


# only needs to be run once for a given cohort (e.g., OneK1K / TOB)
def filter_variants(
    mt_path: str,  # 'mt/v7.mt'
    sc_samples_str: str,
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
    from hail.methods import export_plink

    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)
    logging.info(f'Number of total loci: {mt.count()[0]}')
    logging.info(f'Number of total samples: {mt.count()[1]}')

    # densify
    mt = hl.experimental.densify(mt)

    # add sample filters
    bm_samples = get_bone_marrow_sequencing_groups(mt=mt)
    dup_samples = get_duplicated_samples(mt=mt)
    out_samples = get_non_tob_samples(mt=mt)
    qc_samples = get_low_qc_samples(mt=mt)
    filter_samples = {*bm_samples, *dup_samples, *out_samples, *qc_samples}
    logging.info(f'Total samples to filter: {len(filter_samples)}')
    # use syntax from:
    # https://hail.is/docs/0.2/hail.MatrixTable.html#hail.MatrixTable.filter_cols
    set_to_remove = hl.literal(filter_samples)
    mt = mt.filter_cols(~set_to_remove.contains(mt['s']))
    logging.info(f'Number of samples after filtering: {mt.count()[1]}')
    if mt.count_cols() == 0:
        raise ValueError('No samples left after filtering')

    # subset to relevant samples (samples we have scRNA-seq data for)
    sc_samples = sc_samples_str.split(',')
    logging.info(f'Number of sc samples: {len(sc_samples)}')
    mt = mt.filter_cols(hl.set(sc_samples).contains(mt.s))
    logging.info(f'Number of remaining samples: {mt.count()[1]}')

    # filter out low quality variants and consider biallelic SNPs only
    # (no multi-allelic, no ref-only, no indels)
    mt = mt.filter_rows(
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic (exclude multiallelic)
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
    )

    mt = hl.variant_qc(mt)

    mt.write(output_mt_path, overwrite=True)
    logging.info(f'No QC-passing, biallelic SNPs: {mt.count()[0]}')

    # subset variants for variance ratio estimation
    # minor allele count (MAC) > {vre_n_markers}
    vre_mt = mt.filter_rows(mt.variant_qc.AC[0] > vre_mac_threshold)
    n_ac_vars = vre_mt.count()[0]  # to avoid evaluating this 2X
    logging.info(f'Number of variants post AC filter: {n_ac_vars}')
    if n_ac_vars == 0:
        logging.info('No variants left, exit')
        return
    # since pruning is very costly, subset first a bit
    random.seed(0)
    vre_mt = vre_mt.sample_rows(p=0.01)
    logging.info(f'Initial subset of variants: {vre_mt.count()[0]}')
    # perform LD pruning
    pruned_variant_table = hl.ld_prune(vre_mt.GT, r2=0.2, bp_window_size=500000)
    vre_mt = vre_mt.filter_rows(hl.is_defined(pruned_variant_table[vre_mt.row_key]))
    logging.info(f'Subset of variants after pruning: {vre_mt.count()[0]}')
    # randomly sample {vre_n_markers} variants
    random.seed(0)
    vre_mt = vre_mt.sample_rows((vre_n_markers * 1.1) / vre_mt.count()[0])
    vre_mt = vre_mt.head(vre_n_markers)
    logging.info(f'Subset to {vre_n_markers} variants: {vre_mt.count()[0]}')

    # export to plink common variants only for sparse GRM
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
    # extract CPG IDs from file
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
        sc_samples_str=sc_samples,
        output_mt_path=output_mt_path,
        vre_plink_path=vre_plink_path,
    )


if __name__ == '__main__':
    genotypes_pipeline()  # type: ignore
