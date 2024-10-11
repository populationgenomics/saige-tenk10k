#!/usr/bin/env python3

"""
This script will

- extract sex information
- extract age information
- extract genotype PCs

for the TOB and BioHEART cohorts as part of TenK10K part1

these will be used in the gene expression input preparation
building inputs for the SAIGE-QTL association pipeline.


To run:

analysis-runner \
   --description "get sample covariates" \
   --dataset "bioheart" \
   --access-level "standard" \
   --output-dir "saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/covariates" \
    	python3 get_sample_covariates.py --tenk10k-sampleqc-file-path  'gs://cpg-bioheart-main/large_cohort/tenk10k1-1/sample_qc.ht' \
               --project-names 'tob-wgs,bioheart' --vds-version tenk10k1-1-eur-nbm

"""

from cpg_utils.hail_batch import dataset_path, init_batch, output_path

import click

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
    '--tenk10k-sampleqc-file-path',
    help='this file should contain sample id and inferred sex info for the tenk10k cohort (tob + bioheart)',
)
@click.option('--vds-version', help=' e.g. tenk10k1-0,bioheart1-0 ')
@click.option('--project-names', default='tob-wgs,bioheart')
@click.option('--number-of-sample-perms', default=10)
@click.option('--fill-in-sex', default=False)
@click.option('--fill-in-age', default=False)
@click.command()
def main(
    tenk10k_sampleqc_file_path: str,
    vds_version: str,
    project_names: str,
    number_of_sample_perms: int,
    fill_in_sex: bool,
    fill_in_age: bool,
):
    """
    Get sex, age and genotype PCs for TOB and BioHEART individuals
    """

    init_batch()

    # sex
    # option 1: separate files
    # # tob
    # tob_sample_qc_ht = hl.read_table(tob_sampleqc_file_path)
    # # convert to pandas
    # tob_sample_qc_df = tob_sample_qc_ht.to_pandas()
    # # extract info and reformat to table
    # tob_sample_qc_df['sample_id'] = [str(s) for s in tob_sample_qc_df['s']]
    # # only retain relevant columns
    # tob_sex_df = tob_sample_qc_df[['sample_id', 'sex']]
    # # bioheart
    # bioheart_sample_qc_ht = hl.read_table(bioheart_sampleqc_file_path)
    # # convert to pandas
    # bioheart_sample_qc_df = bioheart_sample_qc_ht.to_pandas()
    # # extract info and reformat to table
    # bioheart_sample_qc_df['sample_id'] = [str(s) for s in bioheart_sample_qc_df['s']]
    # # only retain relevant columns
    # bioheart_sex_df = bioheart_sample_qc_df[['sample_id', 'sex']]
    # # combine_info
    # sex_df = pd.concat([tob_sex_df, bioheart_sex_df], axis=0)

    # # option 2: combined file
    # at the moment this is not up to date, but ideally this would be what we'd run instead
    # sample_qc_ht_path = dataset_path(f'large_cohort/{vds_version}/sample_qc.ht')
    sample_qc_ht = hl.read_table(tenk10k_sampleqc_file_path)
    # convert to pandas
    sample_qc_df = sample_qc_ht.to_pandas()
    # extract info and reformat to table
    sample_qc_df['sample_id'] = [str(s) for s in sample_qc_df['s']]
    # go from sex karyotype (XX, XY) to sex numeric coding (2,1)
    sample_qc_df['sex'] = (
        sample_qc_df['sex_karyotype'].map({"XY": 1, "XX": 2, "X": 0}).astype(int)
    )
    # only retain relevant columns
    sex_df = sample_qc_df[['sample_id', 'sex']]

    # add sex as unknown (0) if missing
    if fill_in_sex:
        sex_df['sex'] = sex_df['sex'].fillna(0)

    # age
    # create a list from dictionary to populate
    age_dict_list = []
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
    # add average age if missing
    if fill_in_age:
        mean_age = age_df['age'].mean()
        age_df['age'].fillna(mean_age)

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
    # drop NAs
    combined_sex_age_pcs.dropna(inplace=True)

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
