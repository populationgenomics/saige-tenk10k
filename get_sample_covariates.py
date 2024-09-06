#!/usr/bin/env python3

"""
This script will

- extract sex information
- extract age information
- extract genotype PCs

for the TOB and BioHEART cohorts as part of TenK10K part1

these will be used in the gene expression input preparation
building inputs for the SAIGE-QTL association pipeline.


To run,

in test:

analysis-runner \
    --description "get sample covariates" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/covariates/" \
    python3 get_sample_covariates.py --tob-sex-file-path 'gs://cpg-bioheart-test/saige-qtl/input_files/mapping_for_anna.csv' \
                --bioheart-sex-file-path 'gs://cpg-bioheart-test-analysis/hoptan-str/somalier/somalier.samples.tsv' \
                --project-names 'tob-wgs,bioheart' --vds-version tenk10k1-0

in main:

analysis-runner \
    --description "get sample covariates" \
    --dataset "bioheart" \
    --access-level "standard" \
    --output-dir "saige-qtl/input_files/covariates/" \
    python3 get_sample_covariates.py --tob-sex-file-path 'gs://cpg-bioheart-test/saige-qtl/input_files/mapping_for_anna.csv' \
                --bioheart-sex-file-path 'gs://cpg-bioheart-main-analysis/qc-stand-alone/somalier/990_samples_somalier.samples.tsv' \
                --project-names 'tob-wgs,bioheart' --vds-version v1-0

"""

from cpg_utils.hail_batch import dataset_path, init_batch, output_path

import click
import sys

import hail as hl
import pandas as pd
from sklearn.utils import shuffle

from metamist.graphql import gql, query

GET_PARTICIPANT_META_QUERY = gql(
    """
    query GetMeta($project_name: String! ) {
        project(name: $project_name) {
            sequencingGroups {
                id
                sample {
                    participant {
                        meta
                    }
                }
            }
        }
    }
    """
)


@click.option(
    '--tob-sex-file-path',
    help='this file should contain sample id and inferred sex info for the tob cohort',
)
@click.option(
    '--bioheart-sex-file-path',
    help='this file should contain sample id and inferred sex info for the bioheart cohort',
)
@click.option('--vds-version', help=' e.g. 1-0 ')
@click.option('--project-names', default='tob-wgs,bioheart')
@click.option('--number-of-sample-perms', default=10)
@click.command()
def main(
    tob_sex_file_path,
    bioheart_sex_file_path,
    vds_version,
    project_names,
    number_of_sample_perms,
):
    """
    Get sex, age and genotype PCs for TOB and BioHEART individuals
    """

    init_batch()

    # sex
    # check if files exist
    try:
        # TOB sex info from Vlad's metadata file
        tob_meta = pd.read_csv(tob_sex_file_path)
        # BioHEART sex info from Hope's Somalier stand alone run
        bioheart_meta = pd.read_csv(bioheart_sex_file_path, sep="\t")
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    # extract sex for TOB
    # remove samples with ambiguous sex inference
    tob_meta = tob_meta[tob_meta['sex_karyotype'].isin(["XX", "XY"])]
    # encode sex as 1,2 instead
    tob_meta['sex'] = tob_meta['sex_karyotype'].replace('XY', '1')
    tob_meta['sex'] = tob_meta['sex'].replace('XX', '2')
    # rename s as sample id to match bioheart file
    tob_meta['sample_id'] = tob_meta['new_CPG_id']
    tob_sex = tob_meta.loc[:, ["sample_id", "sex"]]
    # extract sex for BioHEART
    bioheart_sex = bioheart_meta.loc[:, ["sample_id", "sex"]]
    # combine_info
    sex_df = pd.concat([tob_sex, bioheart_sex], axis=0)

    # age
    # create a list from dictionary to populate
    age_dict_list: list(dict) = []
    # loop over projects (tob-wgs, bioheart)
    for project_name in project_names.split(','):
        query_vars = {'project_name': project_name}
        # run query above, which returns a dict
        meta = query(GET_PARTICIPANT_META_QUERY, variables=query_vars)
        for sg in meta['project']['sequencingGroups']:
            cpg_id = sg['id']
            try:
                age = sg['sample']['participant']['meta']['age']
            except KeyError as e:
                print(f"Key Error: - no {e} available for {cpg_id}")
                age = 'NA'
            age_dict_list.append({'sample_id': cpg_id, 'age': age})
    age_df = pd.DataFrame(age_dict_list)

    # genotype PCs
    pcs_ht_path = dataset_path(f'large_cohort/{vds_version}/ancestry/scores.ht')
    pcs_ht = hl.read_table(pcs_ht_path)
    # convert to pandas
    pcs_df = pcs_ht.to_pandas()
    # extract info and reformat to table
    pcs_df['sample_id'] = [str(s) for s in pcs_df['s']]
    for i in range(len(pcs_df['scores'][0])):
        # extract individual PC values, make numeric after splitting
        pcs_df[f'geno_PC{i+1}'] = [
            float(
                str(score)
                .split(",")[i]
                .replace("[", "")
                .replace("]", "")
                .replace(" ", "")
            )
            for score in pcs_df['scores']
        ]
    # drop unused columns
    pcs_df.drop(columns=['s', 'scores'], inplace=True)

    # index with sample id
    sex_df.index = sex_df['sample_id']
    sex_df.drop(['sample_id'], axis=1, inplace=True)
    age_df.index = age_df['sample_id']
    age_df.drop(['sample_id'], axis=1, inplace=True)
    pcs_df.index = pcs_df['sample_id']
    pcs_df.drop(['sample_id'], axis=1, inplace=True)

    # combine sex, age and geno PC info
    combined_sex_age_pcs = pd.concat([sex_df, age_df, pcs_df], axis=1)

    # add permuted sample ids for calibration analysis
    samples = combined_sex_age_pcs.index
    for i in range(number_of_sample_perms):
        combined_sex_age_pcs[f'sample_perm{i}'] = shuffle(
            samples.values, random_state=i
        )

    # save to file
    sex_age_pcs_out_file = output_path(
        'sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv', 'analysis'
    )
    combined_sex_age_pcs.to_csv(sex_age_pcs_out_file)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
