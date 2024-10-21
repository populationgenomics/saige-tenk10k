#!/usr/bin/env python3

"""
This script will

- open pheno file
- open genotype file
- open conditional files
- create new pheno files with genotype as new column

these files will be used as new inputs for the
SAIGE-QTL conditional analysis pipeline.

To run:

In main:

analysis-runner \
    --description "add variant to pheno files" \
    --dataset "bioheart" \
    --access-level "standard" \
    --output-dir "saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/" \
    python3 make_group_file.py --chromosomes chr21 \
        --pheno-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/pheno_cov_files/ \
        --conditional-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/conditional_files/ \
        --chrom-mt-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/genotypes/vds-tenk10k1-0_qc_pass --max-delay=30

"""

import click
import logging

import hail as hl
import hailtop.batch.job as hb_job

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch


def add_variant_to_pheno_file(
    mt_path: str,
    gene: str,
    chrom: str,
    original_pheno_files_path: str,
    new_pheno_files_path: str,
    genome_reference: str,
):
    """
    Make conditional pheno file
    """
    import random
    import time

    from hail import filter_intervals, parse_locus_interval
    import pandas as pd
    from cpg_utils.hail_batch import init_batch

    init_batch()

    pheno_file = f'{original_pheno_files_path}{chrom}/{gene}.tsv'
    print(f'original pheno file: {pheno_file}')
    gene_df = pd.read_csv(pheno_file, sep='\t')
    num_chrom = gene_df.columns.values[0]
    window_start = gene_df.columns.values[1]
    window_end = gene_df.columns.values[2]
    gene_interval = f'chr{num_chrom}:{window_start}-{window_end}'
    # extract variants within interval
    chrom_mt_filename = f'{mt_path}/{chrom}_rare_variants.mt'
    chrom_mt = hl.read_matrix_table(chrom_mt_filename)
    chrom_mt = filter_intervals(
        chrom_mt,
        [parse_locus_interval(gene_interval, reference_genome=genome_reference)],
    )

    # strip the chr from chromosome, annotate as a new field
    # create a new text field with both alleles
    chrom_mt = chrom_mt.annotate_rows(
        var=hl.delimit(
            [
                chrom_mt.locus.contig.replace('chr', ''),
                hl.str(chrom_mt.locus.position),
                chrom_mt.alleles[0],
                chrom_mt.alleles[1],
            ],
            ':',
        ),
        gene=gene,
    )
    chrom_mt.export(str(group_file).replace('.tsv', '_tmp.tsv'))
    chrom_df = pd.read_csv(str(group_file).replace('.tsv', '_tmp.tsv'), sep='\t')
    chrom_df['anno'] = 'null'
    if gamma != 'none':
        # annos before weights
        chrom_df = chrom_df[['var', 'anno', 'weight:dTSS']]
    vals_df = chrom_df.T
    vals_df['category'] = vals_df.index
    categories_df = pd.DataFrame(categories_data)
    group_vals_df = pd.merge(categories_df, vals_df, on='category')
    with group_file.open('w') as gdf:
        group_vals_df.to_csv(gdf, index=False, header=False, sep=' ')


@click.command()
@click.option('--chromosomes', help=' chr1,chr22 ')
@click.option('--cis-window-files-path')
@click.option('--group-files-path')
@click.option('--chrom-mt-files-path')
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
    chrom_mt_files_path: str,
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
                    mt_path=chrom_mt_files_path,
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