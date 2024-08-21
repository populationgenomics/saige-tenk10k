#!/usr/bin/env python3

"""
This script will

- open genotype file
- open gene cis window files
- create gene rv group files

these files will be used as inputs for the
SAIGE-QTL association pipeline.

To run:

analysis-runner \
    --description "make variant group input files" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/" \
    python3 make_group_file.py --chromosomes chr22 \
        --cis-window-files-path saige-qtl/input_files/ \
        --group-file-path saige-qtl/input_files/group_files/


"""

import click
import logging
import math

import hail as hl
import pandas as pd

import hailtop.batch.job as hb_job
from typing import List

from cpg_utils import to_path
from cpg_utils.config import get_config

from cpg_utils.hail_batch import dataset_path, get_batch, init_batch, output_path


# def distance_to_weight(distance: int, gamma: float = 1e-5):
#     """
#     Define weight for a genetic variant based on
#     the distance of that variant from the gene

#     Following the approach used by the APEX authors
#     doi: https://doi.org/10.1101/2020.12.18.423490
#     """
#     import math

#     weight = math.exp(-gamma * abs(distance))
#     return weight


def build_group_file_single_gene(gene: str, out_path: str, variants, weights):
    """
    Build group file for SAIGE-QTL

    First row: genetic variants to test
    Second row: annotations (none here)
    Third row: weights (dTSS here)
    """
    group_df = pd.DataFrame(
        {'gene': [gene, gene, gene], 'category': ['var', 'anno', 'weight:dTSS']}
    )
    vals_df = pd.DataFrame({'var': variants, 'anno': 'null', 'weight:dTSS': weights}).T
    vals_df['category'] = vals_df.index
    # combine
    group_vals_df = pd.merge(group_df, vals_df, on='category')
    print(group_vals_df.head())
    print(f'out path: {out_path}')
    print(f'output_path(out_path): {output_path(out_path)}')
    with to_path(output_path(out_path)).open('w') as gdf:
        group_vals_df.to_csv(gdf, index=False, header=False)


@click.command()
@click.option('--chromosomes', help=' chr1,chr22 ')
@click.option('--cis-window-files-path')
@click.option('--group-file-path')
@click.option('--cis-window', default=100000)
@click.option('--gamma', default=1e-5)
@click.option(
    '--concurrent-job-cap',
    type=int,
    default=100,
    help=(
        'To avoid resource starvation, set this concurrency to limit '
        'horizontal scale. Higher numbers have a better walltime, but '
        'risk jobs that are stuck (which are expensive)'
    ),
)
def main(
    chromosomes: str,
    cis_window_files_path: str,
    group_file_path: str,
    cis_window: int,
    gamma: float,
    concurrent_job_cap: int,
):
    """
    Run expression processing pipeline
    """
    # set this up with the default (scanpy) python image
    get_batch(
        default_python_image=get_config()['images']['scanpy'],
        name='prepare all gene files',
    )
    all_jobs: List[hb_job.Job] = []

    def manage_concurrency(new_job: hb_job.Job):
        """
        Manage concurrency, so that there is a cap on simultaneous jobs
        Args:
            new_job (hb_job.Job): a new job to add to the stack
        """
        if len(all_jobs) > concurrent_job_cap:
            new_job.depends_on(all_jobs[-concurrent_job_cap])
        all_jobs.append(new_job)

    init_batch()

    # loop over chromosomes
    for chrom in chromosomes.split(','):
        print(f'chrom: {chrom}')

        # load rare variant vcf file for specific chromosome
        vcf_path = dataset_path(
            f'saige-qtl/bioheart/input_files/genotypes/vds-bioheart1-0/{chrom}_rare_variants.vcf.bgz'
        )
        ds = hl.import_vcf(vcf_path, reference_genome='GRCh37')

        # print(cis_window_files_path)
        # do a glob, then pull out all file names as Strings
        files = [
            str(file)
            for file in to_path(cis_window_files_path).glob(
                f'cis_window_files/{chrom}/*bp.tsv'
            )
        ]
        # test only
        files = files[0:10]
        logging.info(f'I found these files: {", ".join(files)}')

        genes = [
            f.replace(f'_{cis_window}bp.tsv', '').replace(
                f'{cis_window_files_path}cis_window_files/{chrom}/', ''
            )
            for f in files
        ]
        logging.info(f'I found these genes: {", ".join(genes)}')

        for gene in genes:
            print(f'gene: {gene}')
            # get gene cis window info
            gene_file = to_path(
                f'{cis_window_files_path}cis_window_files/{chrom}/{gene}_{cis_window}bp.tsv'
            )
            print(f'gene file: {gene_file}')
            gene_df = pd.read_csv(gene_file, sep='\t')
            num_chrom = gene_df.columns.values[0]
            window_start = gene_df.columns.values[1]
            window_end = gene_df.columns.values[2]
            gene_interval = f'{num_chrom}:{window_start}-{window_end}'
            # extract variants within interval
            ds_result = hl.filter_intervals(
                ds, [hl.parse_locus_interval(gene_interval, reference_genome='GRCh37')]
            )
            variants = [loc.position for loc in ds_result.locus.collect()]
            gene_tss = int(window_start) + cis_window
            distances = [int(var) - gene_tss for var in variants]
            # get weight for genetic variants based on
            # the distance of that variant from the gene
            # Following the approach used by the APEX authors
            # doi: https://doi.org/10.1101/2020.12.18.423490
            weights = [math.exp(-gamma * abs(d)) for d in distances]
            print(f'group file: {group_file_path}')
            group_file_job = get_batch().new_python_job(name=f'group file: {gene}')
            group_file_job.call(
                build_group_file_single_gene,
                gene,
                str(group_file_path),
                variants,
                weights,
            )
            manage_concurrency(group_file_job)

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
