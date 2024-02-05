#!/usr/bin/env python3

"""
This script will

- extract sex information
- extract age information

for the TOB and BioHEART cohorts as part of TenK10K part1

these will be used in the gene expression input preparation
building inputs for the SAIGE-QTL association pipeline.


To run:

analysis-runner \
    --description "get sample covariates" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/covariates/" \
    python3 get_sample_covariates.py --tob-sex-file-path 'gs://cpg-tob-wgs-test-analysis/joint-calling/v7/meta.tsv' \
                --bioheart-sex-file-path 'gs://cpg-bioheart-test-analysis/hoptan-str/somalier/somalier.samples.tsv'

main files:
'gs://cpg-tob-wgs-main-analysis/joint-calling/v7/meta.tsv'
'gs://cpg-bioheart-main-analysis/qc-stand-alone/somalier/990_samples_somalier.samples.tsv
"""

# from cpg_utils import to_path
# from cpg_utils.config import get_config
# from cpg_utils.hail_batch import dataset_path, get_batch, init_batch, output_path
from cpg_utils.hail_batch import output_path

import click
import pandas as pd


@click.option('--tob-sex-file-path')
@click.option('--bioheart-sex-file-path')
@click.command()
def main(tob_sex_file_path, bioheart_sex_file_path):
    """
    Get sex and age info for TOB and BioHEART individuals
    """
    # TOB sex info from Vlad's metadata file
    tob_meta = pd.read_csv(tob_sex_file_path, sep="\t")
    # remove non-TOB samples
    tob_meta = tob_meta[
        ~tob_meta['s'].isin(["NA12878", "NA12891", "NA12892", "syndip"])
    ]
    print(tob_meta['sex_karyotype'].unique())
    # remove samples with ambiguous sex inference
    tob_meta = tob_meta[tob_meta['sex_karyotype'].isin(["XX", "XY"])]
    # encode sex as 1,2 instead
    tob_meta['sex'] = tob_meta['sex_karyotype'].replace('XY', '1')
    tob_meta['sex'] = tob_meta['sex_karyotype'].replace('XX', '2')
    # rename s as sample id to match bioheart file
    tob_meta['sample_id'] = tob_meta['s']
    tob_sex = tob_meta.loc[:, ["sample_id", "sex"]]
    print(tob_sex['sex'].unique())
    # BioHEART sex info from Hope's Somalier stand alone run
    bioheart_meta = pd.read_csv(bioheart_sex_file_path, sep="\t")
    bioheart_sex = bioheart_meta.loc[:, ["sample_id", "sex"]]
    print(bioheart_sex['sex'].unique())
    # combine_info
    combined_sex = pd.concat([tob_sex, bioheart_sex], axis=0)
    sex_out_file = output_path('sex_tob_bioheart.csv')
    combined_sex.to_csv(sex_out_file)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
