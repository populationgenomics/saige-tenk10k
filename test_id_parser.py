#!/usr/bin/env python3
# pylint: disable=import-error,line-too-long


"""

Replaces CPG IDs in pheno cov files with test CPG ids for testing purposes

analysis-runner --dataset "bioheart" \
    --description "Replace CPG IDs in pheno cov files with test CPG ids for testing purposes" \
    --access-level "test" \
    --output-dir "str/saige-qtl/input_files" \
    test_id_parser.py --pheno-cov-path=gs://cpg-bioheart-test/saige-qtl/input_files/pheno_cov_files
    --cell-type=B_naive --chrom=chr22

"""
import pandas as pd
import click
from cpg_utils import to_path

from cpg_utils.hail_batch import get_batch


def id_parser(pheno_cov_path: str, gcs_output_path: str):
    """
    Replaces CPG IDs in pheno cov files with test CPG ids for testing purposes
    """
    data = pd.read_csv(pheno_cov_path, sep='\t')
    individuals_to_filter = [
        'CPG247841',
        'CPG247361',
        'CPG247775',
        'CPG247643',
        'CPG247502',
    ]
    filtered_df = data[data['individual'].isin(individuals_to_filter)]

    # Mapping to rename individuals
    mapping = {
        'CPG247841': 'CPG305235',
        'CPG247361': 'CPG305326',
        'CPG247775': 'CPG305359',
        'CPG247643': 'CPG305342',
        'CPG247502': 'CPG305334',
    }

    # Rename individuals based on the mapping
    filtered_df['individual'] = filtered_df['individual'].replace(mapping)

    filtered_df.to_csv(gcs_output_path, sep='\t', index=False)


@click.option(
    '--pheno-cov-path',
    help='GCS file path to pheno cov files',
    type=str,
)
@click.option(
    '--cell-type',
    help='Cell type to filter',
    type=str,
)
@click.option(
    '--chrom',
    help='Chromosome to filter',
    type=str,
)
@click.command()
def main(pheno_cov_path: str, cell_type: str, chrom: str):
    files = list(
        to_path(
            f'{pheno_cov_path}/{cell_type}/{chrom}',
        ).glob('*.tsv')
    )
    for file in files:
        file_path = str(file)
        output_path = f'{pheno_cov_path}/remapped_ids/{cell_type}/{chrom}/{file_path}'
        id_parser(file_path, output_path)



if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
