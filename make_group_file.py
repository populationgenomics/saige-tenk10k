#!/usr/bin/env python3

"""
This script will

- open genotype file
- open gene cis window files
- create gene rv group files

these files will be used as inputs for the
SAIGE-QTL association pipeline.

To run:

In main:

analysis-runner \
    --description "make variant group input files" \
    --dataset "tenk10k" \
    --access-level "full" \
    --output-dir "saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/" \
    python3 make_group_file.py --chromosomes chr22 \
        --cis-window-files-path=gs://cpg-tenk10k-main/saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/cis_window_files/ \
        --group-files-path=gs://cpg-tenk10k-main/saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/group_files_mt_n1925/annotated/ \
        --chrom-mt-files-path=gs://cpg-tenk10k-main/saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/genotypes/vds-tenk10k1-0_qc_pass --max-delay=30 \
        --rv-annotation-files-path=gs://cpg-tenk10k-test/saige-qtl/rv_annot/

In test:

analysis-runner \
   --description "make variant group input files" \
   --dataset "bioheart" \
   --access-level "test" \
   --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/input_files/" \
   python3 make_group_file.py --chromosomes chr21 \
       --cis-window-files-path gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/input_files/cis_window_files/ \
       --group-files-path gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/input_files/group_files_mt/ \
       --chrom-mt-files-path gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/input_files/genotypes/vds-tenk10k1-0_subset --max-delay=10

"""

import click
import logging

import hail as hl
import hailtop.batch.job as hb_job

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch


def make_group_file(
    mt_path: str,
    gene: str,
    chrom: str,
    cis_window_files_path: str,
    group_file,
    rare_variant_annotation_files_path: str,
    cis_window: int,
    genome_reference: str,
    gamma: str,
    max_delay: int,
):
    """
    Make group file
    """
    import random
    import time

    from hail import filter_intervals, parse_locus_interval
    import numpy as np
    import pandas as pd
    from cpg_utils.hail_batch import init_batch

    time.sleep(random.randint(0, max_delay))

    init_batch()

    gene_file = f'{cis_window_files_path}{chrom}/{gene}_{cis_window}bp.tsv'
    print(f'gene file: {gene_file}')
    gene_df = pd.read_csv(gene_file, sep='\t')
    num_chrom = gene_df.columns.values[0]
    window_start = round(float(gene_df.columns.values[1]))
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

    # add distance-based weights if required

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
            **{'weight:dTSS': hl.exp(-float(gamma) * chrom_mt.distance)}
        )

        # re-key by variant string instead of locus/alleles
        # select the columns we want (dropping the rest)
        # keep only the rows
        chrom_mt = (
            chrom_mt.key_rows_by(chrom_mt.var).select_rows('gene', 'weight:dTSS').rows()
        )
        categories_data = {
            'gene': [gene, gene, gene],
            'category': ['var', 'anno', 'weight:dTSS'],
        }

    else:
        # we don't annotate weights/distances
        chrom_mt = chrom_mt.key_rows_by(chrom_mt.var).select_rows('gene').rows()
        categories_data = {'gene': [gene, gene, gene], 'category': ['var', 'anno']}

    chrom_mt.export(str(group_file).replace('.tsv', '_tmp.tsv'))
    chrom_df = pd.read_csv(str(group_file).replace('.tsv', '_tmp.tsv'), sep='\t')

    # add annotations if required

    if rare_variant_annotation_files_path != 'none':
        annot_file = (
            f'{rare_variant_annotation_files_path}{chrom}_rare_variant_annotations.csv'
        )
        annot_df = pd.read_csv(annot_file)
        # list the functional columns you care about (only keep those that actually exist)
        requested_functional = [
            "is_promoter",
            "is_enhancer",
            "is_ctcf",
            "is_dnase",
            "is_tx",
            "is_atac",
            "is_open",
            "is_tf",
            "is_promoter_ai_gte0.2",
            "is_splice_ai_gte0.2",
        ]
        functional_cols = [c for c in requested_functional if c in annot_df.columns]

        # find all current "is_" columns (may or may not include is_functional yet)
        all_is_cols = [c for c in annot_df.columns if c.startswith("is_")]

        # convert all existing is_ columns EXCEPT is_functional to bool (fill NAs -> False)
        cols_to_bool = [c for c in all_is_cols if c != "is_functional"]
        annot_df[cols_to_bool] = annot_df[cols_to_bool].fillna(False).astype(bool)

        # compute is_functional from the cleaned functional_cols (safe because they are bool now)
        if functional_cols:
            annot_df["is_functional"] = annot_df[functional_cols].any(axis=1)
        else:
            annot_df["is_functional"] = False

        # rebuild the list of is_ columns (this will include is_functional)
        anno_cols = [c for c in annot_df.columns if c.startswith("is_")]

        # build the anno string from the boolean is_ columns (will include "functional" if True)
        annot_df["anno"] = annot_df[anno_cols].apply(
            lambda row: ",".join([c[3:] for c in anno_cols if row[c]]), axis=1
        )
        annot_df_sel = annot_df[['POS', 'anno']]
        # Merge on POS
        combined_df = chrom_df.merge(annot_df_sel, on='POS', how='left')
    # if no annotations are provided
    else:
        combined_df = chrom_df
        combined_df['anno'] = 'null'

    if gamma != 'none':
        # annos before weights
        combined_df = combined_df[['var', 'anno', 'weight:dTSS']]
    vals_df = combined_df.T
    vals_df['category'] = vals_df.index
    categories_df = pd.DataFrame(categories_data)
    group_vals_df = pd.merge(categories_df, vals_df, on='category')
    # drop any variants with no annotations (after making all nas)
    group_vals_df = group_vals_df.replace(r'^\s*$', np.nan, regex=True)
    group_vals_df_annos = group_vals_df.dropna(axis=1)
    with group_file.open('w') as gdf:
        group_vals_df_annos.to_csv(gdf, index=False, header=False, sep=' ')


@click.command()
@click.option('--chromosomes', help=' chr1,chr22 ')
@click.option('--cis-window-files-path')
@click.option('--group-files-path')
@click.option('--chrom-mt-files-path')
@click.option('--rare-variant-annotation-files-path')
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
    rare_variant_annotation_files_path: str,
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
                    rare_variant_annotation_files_path=rare_variant_annotation_files_path,
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
