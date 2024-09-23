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
        --cis-window-files-path gs://cpg-bioheart-test/saige-qtl/input_files/cis_window_files/ \
        --group-files-path gs://cpg-bioheart-test/saige-qtl/input_files/group_files/ \
        --vcf-path gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/input_files/genotypes/v3/vds-bioheart1-0 \
        --ngenes-to-test 5


"""

import click
import logging

# import math

import hail as hl

# import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch


def make_group_file(
    # ds: hl.vds,
    vcf_path_chrom: str,
    gene: str,
    chrom: str,
    cis_window_files_path: str,
    # group_files_path: str,
    group_file,
    cis_window: int,
    genome_reference: str,
    gamma: str,
):
    """
    Make group file
    """
    import math

    # import hail as hl
    from hail import filter_intervals, import_vcf, parse_locus_interval
    import pandas as pd

    # from cpg_utils import to_path

    gene_file = f'{cis_window_files_path}{chrom}/{gene}_{cis_window}bp.tsv'
    print(f'gene file: {gene_file}')
    gene_df = pd.read_csv(gene_file, sep='\t')
    num_chrom = gene_df.columns.values[0]
    window_start = gene_df.columns.values[1]
    window_end = gene_df.columns.values[2]
    gene_interval = f'{num_chrom}:{window_start}-{window_end}'
    # extract variants within interval
    # ds_result = hl.filter_intervals(
    #     ds,
    #     [hl.parse_locus_interval(gene_interval, reference_genome=genome_reference)],
    # )
    ds = import_vcf(vcf_path_chrom, reference_genome=genome_reference)
    ds_result = filter_intervals(
        ds,
        [parse_locus_interval(gene_interval, reference_genome=genome_reference)],
    )
    variants_chrom_pos = [
        f'{loc.contig}:{loc.position}' for loc in ds_result.locus.collect()
    ]
    variants_alleles = [
        f'{allele[0]}:{allele[1]}' for allele in ds_result.alleles.collect()
    ]
    variants = [
        f'{variants_chrom_pos[i]}:{variants_alleles[i]}'
        for i in range(len(variants_chrom_pos))
    ]

    if gamma != 'none':
        gene_tss = int(window_start) + cis_window
        distances = [int(var.split(":")[1]) - gene_tss for var in variants]
        # get weight for genetic variants based on
        # the distance of that variant from the gene
        # Following the approach used by the APEX authors
        # doi: https://doi.org/10.1101/2020.12.18.423490
        weights = [math.exp(-float(gamma) * abs(d)) for d in distances]
        group_df = pd.DataFrame(
            {
                'gene': [gene, gene, gene],
                'category': ['var', 'anno', 'weight:dTSS'],
            }
        )
        vals_df = pd.DataFrame(
            {'var': variants, 'anno': 'null', 'weight:dTSS': weights}
        ).T
        # group_file = f'{group_files_path}{chrom}/{gene}_{cis_window}bp_dTSS_weights.tsv'
    else:
        group_df = pd.DataFrame({'gene': [gene, gene], 'category': ['var', 'anno']})
        vals_df = pd.DataFrame({'var': variants, 'anno': 'null'}).T
        # group_file = f'{group_files_path}{chrom}/{gene}_{cis_window}bp_no_weights.tsv'
    vals_df['category'] = vals_df.index
    # combine
    group_vals_df = pd.merge(group_df, vals_df, on='category')

    # with to_path(group_file).open('w') as gdf:
    #     group_vals_df.to_csv(gdf, index=False, header=False, sep=' ')
    with group_file.open('w') as gdf:
        group_vals_df.to_csv(gdf, index=False, header=False, sep=' ')


@click.command()
@click.option('--chromosomes', help=' chr1,chr22 ')
@click.option('--cis-window-files-path')
@click.option('--group-files-path')
@click.option('--vcf-path')
@click.option('--cis-window', default=100000)
@click.option('--gamma', default='1e-5')
@click.option('--ngenes-to-test', default='all')
@click.option('--genome-reference', default='GRCh37')
def main(
    chromosomes: str,
    cis_window_files_path: str,
    group_files_path: str,
    vcf_path: str,
    cis_window: int,
    gamma: str,
    ngenes_to_test: str,
    genome_reference: str,
):
    """
    Make group file for rare variant pipeline
    """

    init_batch()

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

        # load rare variant vcf file for specific chromosome
        vcf_path_chrom = f'{vcf_path}/{chrom}_rare_variants.vcf.bgz'
        # ds = hl.import_vcf(vcf_path_chrom, reference_genome=genome_reference)

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
            gene_group_job = get_batch().new_python_job(
                name=f'gene make group file: {gene}'
            )
            # gene_group_job.cpu(group_job_cpu)
            gene_group_job.call(
                make_group_file,
                # ds=ds,
                vcf_path_chrom=vcf_path_chrom,
                gene=gene,
                chrom=chrom,
                cis_window_files_path=cis_window_files_path,
                # group_files_path=group_files_path,
                group_file=to_path(group_file),
                cis_window=cis_window,
                genome_reference=genome_reference,
                gamma=gamma,
            )
            # manage_concurrency(gene_cis_job)
            logging.info(f'make group file job for {gene} scheduled')
            # # get gene cis window info
            # gene_file = f'{cis_window_files_path}{chrom}/{gene}_{cis_window}bp.tsv'
            # print(f'gene file: {gene_file}')
            # gene_df = pd.read_csv(gene_file, sep='\t')
            # num_chrom = gene_df.columns.values[0]
            # window_start = gene_df.columns.values[1]
            # window_end = gene_df.columns.values[2]
            # gene_interval = f'{num_chrom}:{window_start}-{window_end}'
            # # extract variants within interval
            # ds_result = hl.filter_intervals(
            #     ds,
            #     [
            #         hl.parse_locus_interval(
            #             gene_interval, reference_genome=genome_reference
            #         )
            #     ],
            # )
            # variants_chrom_pos = [
            #     f'{loc.contig}:{loc.position}' for loc in ds_result.locus.collect()
            # ]
            # variants_alleles = [
            #     f'{allele[0]}:{allele[1]}' for allele in ds_result.alleles.collect()
            # ]
            # variants = [
            #     f'{variants_chrom_pos[i]}:{variants_alleles[i]}'
            #     for i in range(len(variants_chrom_pos))
            # ]

            # if gamma != 'none':
            #     gene_tss = int(window_start) + cis_window
            #     distances = [int(var.split(":")[1]) - gene_tss for var in variants]
            #     # get weight for genetic variants based on
            #     # the distance of that variant from the gene
            #     # Following the approach used by the APEX authors
            #     # doi: https://doi.org/10.1101/2020.12.18.423490
            #     weights = [math.exp(-float(gamma) * abs(d)) for d in distances]
            #     group_df = pd.DataFrame(
            #         {
            #             'gene': [gene, gene, gene],
            #             'category': ['var', 'anno', 'weight:dTSS'],
            #         }
            #     )
            #     vals_df = pd.DataFrame(
            #         {'var': variants, 'anno': 'null', 'weight:dTSS': weights}
            #     ).T
            #     group_file = (
            #         f'{group_files_path}{chrom}/{gene}_{cis_window}bp_dTSS_weights.tsv'
            #     )
            # else:
            #     group_df = pd.DataFrame(
            #         {'gene': [gene, gene], 'category': ['var', 'anno']}
            #     )
            #     vals_df = pd.DataFrame({'var': variants, 'anno': 'null'}).T
            #     group_file = (
            #         f'{group_files_path}{chrom}/{gene}_{cis_window}bp_no_weights.tsv'
            #     )
            # vals_df['category'] = vals_df.index
            # # combine
            # group_vals_df = pd.merge(group_df, vals_df, on='category')

            # with to_path(group_file).open('w') as gdf:
            #     group_vals_df.to_csv(gdf, index=False, header=False, sep=' ')
    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
