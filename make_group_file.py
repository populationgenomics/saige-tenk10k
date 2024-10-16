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
   --access-level "standard" \
   --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920" \
   python3 make_group_file.py --chromosomes chr2 \
       --cis-window-files-path gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/cis_window_files/ \
       --group-files-path gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/group_files/ \
       --vcf-path gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/genotypes/vds-tenk10k1-0

"""

import click
import logging

import hail as hl
import hailtop.batch.job as hb_job

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch


def make_group_file(
    vds_path: str,
    gene: str,
    chrom: str,
    cis_window_files_path: str,
    group_file,
    cis_window: int,
    genome_reference: str,
    gamma: str,
    max_delay: int,
):
    """
    Make group file
    """
    import math
    import random
    import time

    from hail import filter_intervals, parse_locus_interval
    import pandas as pd
    from cpg_utils.hail_batch import init_batch

    time.sleep(random.randint(0, max_delay))

    init_batch()

    gene_file = f'{cis_window_files_path}{chrom}/{gene}_{cis_window}bp.tsv'
    print(f'gene file: {gene_file}')
    gene_df = pd.read_csv(gene_file, sep='\t')
    num_chrom = gene_df.columns.values[0]
    window_start = gene_df.columns.values[1]
    window_end = gene_df.columns.values[2]
    gene_interval = f'chr{num_chrom}:{window_start}-{window_end}'
    # extract variants within interval - ideally drop this read & densify (see #193)
    vds = hl.vds.read_vds(vds_path)
    chrom_vds = hl.vds.filter_chromosomes(vds, keep=chrom)
    chrom_mt = hl.vds.to_dense_mt(chrom_vds)
    chrom_mt = filter_intervals(
        chrom_mt,
        [parse_locus_interval(gene_interval, reference_genome=genome_reference)],
    )

    # strip the chr from chromosome, annotate as a new field
    # create a new text field with both alleles
    chrom_mt = chrom_mt.annotate_rows(
        variant_string=hl.delimit(
            [
                chrom_mt.chrom.replace('chr', ''),
                hl.str(chrom_mt.locus.position),
                *chrom_mt.alleles,
            ],
            ':',
        ),
        gene=gene,
    )

    if gamma != 'none':

        gene_tss = int(window_start) + cis_window

        # annotate distances
        chrom_mt = chrom_mt.annotate_rows(
            distance=hl.abs(chrom_mt.locus.position - gene_tss)
        )

        # get weight for genetic variants based on
        # the distance of that variant from the gene
        # Following the approach used by the APEX authors
        # doi: https://doi.org/10.1101/2020.12.18.423490
        chrom_mt = chrom_mt.annotate_rows(
            **{'weight:dTSS': math.exp(-float(gamma) * chrom_mt.distance)}
        )

        # re-key by variant string instead of locus/alleles
        # select the columns we want (dropping the rest)
        # keep only the rows
        chrom_mt = (
            chrom_mt.key_rows_by(chrom_mt)
            .select_rows('variant_string', 'gene', 'distance', 'weight:dTSS')
            .rows()
        )

    else:
        # we don't annotate weights/distances
        chrom_mt = (
            chrom_mt.key_rows_by(chrom_mt).select_rows('variant_string', 'gene').rows()
        )

    chrom_mt.export(group_file)


@click.command()
@click.option('--chromosomes', help=' chr1,chr22 ')
@click.option('--cis-window-files-path')
@click.option('--group-files-path')
@click.option('--vds-path')
@click.option('--cis-window', default=100000)
@click.option('--gamma', default='1e-5')
@click.option('--ngenes-to-test', default='all')
@click.option('--genome-reference', default='GRCh38')
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
@click.option(
    '--max-delay',
    type=int,
    default=3000,
    help='delay starting the jobs as they all access the same VDS which causes Hail issues',
)
@click.option(
    '--gene-group-storage',
    default='8G',
)
@click.option(
    '--gene-group-memory',
    default='8G',
)
def main(
    chromosomes: str,
    cis_window_files_path: str,
    group_files_path: str,
    vds_path,
    cis_window: int,
    gamma: str,
    ngenes_to_test: str,
    genome_reference: str,
    concurrent_job_cap: int,
    max_delay: int,
    gene_group_storage: str,
    gene_group_memory: str,
):
    """
    Make group file for rare variant pipeline
    """

    init_batch()

    all_jobs: list[hb_job.Job] = []

    def manage_concurrency(new_job: hb_job.Job):
        """
        Manage concurrency, so that there is a cap on simultaneous jobs
        Args:
            new_job (hb_job.Job): a new job to add to the stack
        """
        if len(all_jobs) > concurrent_job_cap:
            new_job.depends_on(all_jobs[-concurrent_job_cap])
        all_jobs.append(new_job)

    # loop over chromosomes
    for chrom in chromosomes.split(','):
        print(f'chrom: {chrom}')

        # do a glob, then pull out all file names as Strings
        files = [
            str(file)
            for file in to_path(cis_window_files_path).glob(f'{chrom}/*bp.tsv')
        ]
        # if specified, only test ngenes genes
        if ngenes_to_test != 'all':
            files = files[0 : int(ngenes_to_test)]
        logging.info(f'I found these files: {", ".join(files)}')

        genes = [
            f.replace(f'_{cis_window}bp.tsv', '').replace(
                f'{cis_window_files_path}{chrom}/', ''
            )
            for f in files
        ]
        logging.info(f'I found these genes: {", ".join(genes)}')

        for gene in genes:
            print(f'gene: {gene}')
            if gamma != 'none':
                group_file = (
                    f'{group_files_path}{chrom}/{gene}_{cis_window}bp_dTSS_weights.tsv'
                )
            else:
                group_file = (
                    f'{group_files_path}{chrom}/{gene}_{cis_window}bp_no_weights.tsv'
                )
            if not to_path(group_file).exists():
                gene_group_job = get_batch().new_python_job(
                    name=f'gene make group file: {gene}'
                )
                gene_group_job.storage(gene_group_storage)
                gene_group_job.memory(gene_group_memory)
                gene_group_job.call(
                    make_group_file,
                    vds_path=vds_path,
                    gene=gene,
                    chrom=chrom,
                    cis_window_files_path=cis_window_files_path,
                    group_file=to_path(group_file),
                    cis_window=cis_window,
                    genome_reference=genome_reference,
                    gamma=gamma,
                    max_delay=max_delay,
                )
                manage_concurrency(gene_group_job)
                logging.info(f'make group file job for {gene} scheduled')

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
