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
    python3 make_group_file.py --chromosomes chr1


"""

import click
import logging

import hailtop.batch.job as hb_job

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch

@click.command()
@click.option('--cis-window-files-path')
@click.option('--cis-window', default=100000)
# @click.option(
#     '--concurrent-job-cap',
#     type=int,
#     default=100,
#     help=(
#         'To avoid resource starvation, set this concurrency to limit '
#         'horizontal scale. Higher numbers have a better walltime, but '
#         'risk jobs that are stuck (which are expensive)'
#     ),
# )
def main(
    cis_window_files_path: str,
    cis_window: int,
    # concurrent_job_cap: int,
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

    # def manage_concurrency(new_job: hb_job.Job):
    #     """
    #     Manage concurrency, so that there is a cap on simultaneous jobs
    #     Args:
    #         new_job (hb_job.Job): a new job to add to the stack
    #     """
    #     if len(all_jobs) > concurrent_job_cap:
    #         new_job.depends_on(all_jobs[-concurrent_job_cap])
    #     all_jobs.append(new_job)

    init_batch()

    # do a glob, then pull out all file names as Strings
    files = [
        file.name
        for file in to_path(cis_window_files_path).glob(
            'cis_window_files/*/*bp.csv'
        )
    ]
    logging.info(f'I found these files: {", ".join(files)}')

    genes = [f.replace(f'_{cis_window}bp.tsv', '') for f in files]
    logging.info(f'I found these genes: {", ".join(genes)}')


    # for gene in genes:


    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter