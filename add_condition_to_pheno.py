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
    python3 add_condition_to_pheno.py \
        --pheno-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/pheno_cov_files/ \
        --conditional-files-path=gs://cpg-bioheart-main-analysis/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/conditional_files/ \
        --chrom-mt-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/input_files/genotypes/vds-tenk10k1-0_qc_pass --max-delay=30

In test:

analysis-runner \
   --description "add variant to pheno files" \
   --dataset "bioheart" \
   --access-level "test" \
   --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/input_files/" \
   python3 add_condition_to_pheno.py \
       --pheno-files-path gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/input_files/pheno_cov_files/ \
       --conditional-files-path gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100// \
       --chrom-mt-files-path gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/input_files/genotypes/vds-tenk10k1-0_subset --max-delay=10


"""

import click
import logging
import pandas as pd

import hail as hl
import hailtop.batch.job as hb_job

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch


def add_variant_to_pheno_file(
    mt_path: str,
    gene: str,
    chrom: str,
    celltype: str,
    original_pheno_files_path: str,
    new_pheno_files_path,
    conditional_files_path: str,
):
    """
    Add variant as column to pheno cov file
    """
    import pandas as pd
    from cpg_utils.hail_batch import init_batch

    init_batch()

    # open original pheno cov file
    pheno_file = (
        f'{original_pheno_files_path}{celltype}/{chrom}/{gene}_{celltype}_pheno_cov.tsv'
    )
    print(f'original pheno file: {pheno_file}')
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    # open conditional file to extract variant(s)
    conditional_file = f'{conditional_files_path}{celltype}_conditional_file.tsv'
    condition_df = pd.read_csv(conditional_file, sep='\t')
    variant = condition_df[condition_df['gene'] == gene]['variants_to_condition_on']
    # extract variant(s) from chrom mt
    chrom_mt_filename = f'{mt_path}/{chrom}_common_variants.mt'
    chrom_mt = hl.read_matrix_table(chrom_mt_filename)
    # extract genotypes for relevant variant(s)
    chrom_mt_filtered = chrom_mt.filter_rows(
        (chrom_mt.locus.position == int(variant.split(':')[1]))
        & (chrom_mt.alleles[1] == variant.split(':')[3])
    )
    # steal all the entries as a Table, dropping everything except chr, pos, alleles, Genotypes
    genos = chrom_mt_filtered.select_entries('GT').select_rows().entries()
    # create an integer representation of the genotypes
    genos = genos.annotate(
        GT=hl.case()
        .when(genos.GT.is_hom_var(), 2)
        .when(genos.GT.is_het(), 1)
        .default(0)
    )
    # export
    genos.export('table.tsv', delimiter='\t')
    variant_underscores = variant.replace(":", "_")
    new_pheno_file = pheno_file.replace('.tsv', f'_{variant_underscores}.tsv')
    genos.export(str(new_pheno_file).replace('.tsv', '_tmp.tsv'), delimiter='\t')
    geno_df = pd.read_csv(str(new_pheno_file).replace('.tsv', '_tmp.tsv'), sep='\t')
    # rename useful columns and drop the rest
    geno_df = geno_df.rename(columns={"s": "individual", "GT": variant})
    geno_df = geno_df.drop(['locus', 'alleles'], axis=1)

    # add as column to new df (merging on pheno file)
    new_pheno_df = pd.merge(geno_df, pheno_df, on='individual', how='right')

    with new_pheno_files_path.open('w') as npf:
        new_pheno_df.to_csv(npf, index=False, header=False, sep=' ')


@click.command()
@click.option('--pheno-files-path')
@click.option('--condition-pheno-files-path')
@click.option('--conditional-files-path')
@click.option('--chrom-mt-files-path')
@click.option('--ngenes-to-test', default='all')
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
    '--gene-condition-storage',
    default='8G',
)
@click.option(
    '--gene-condition-memory',
    default='8G',
)
def main(
    pheno_files_path: str,
    condition_pheno_files_path: str,
    conditional_files_path: str,
    chrom_mt_files_path: str,
    ngenes_to_test: str,
    concurrent_job_cap: int,
    gene_condition_storage: str,
    gene_condition_memory: str,
):
    """
    Make conditional pheno cov file for conditional analysis
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

    # extract conditional files using glob
    files = [
        str(file)
        for file in to_path(conditional_files_path).glob('*_conditional_file.tsv')
    ]
    # determine what cell type we have conditional files for
    celltypes = [
        file.replace(conditional_files_path, '').replace('_conditional_file.tsv', '')
        for file in files
    ]

    # loop over celltypes
    for celltype in celltypes:
        print(f'celltype: {celltype}')
        # pheno cov file path
        pheno_files_path_ct = f'{pheno_files_path}{celltype}/'
        # open cell type specific conditional file
        conditional_files_path_ct_file = (
            f'{conditional_files_path}{celltype}_conditional_file.tsv'
        )
        conditional_df = pd.read_csv(conditional_files_path_ct_file, sep='\t')
        # extract unique chromosomes
        conditional_df['chr'] = [
            f"chr{variant.split(':')[0]}" for variant in conditional_df['top_MarkerID']
        ]
        # loop over chromosomes
        for chrom in conditional_df['chr'].unique():
            print(f'chrom: {chrom}')

            conditional_df_chr = conditional_df[conditional_df['chr'] == chrom]
            genes = conditional_df_chr['gene']
            logging.info(f'genes to test: {", ".join(genes)}')

            # if specified, only test ngenes genes
            if ngenes_to_test != 'all':
                genes = genes[0 : int(ngenes_to_test)]
            logging.info(f'I found these files: {", ".join(files)}')

            for gene in genes:
                print(f'gene: {gene}')
                pheno_original_file = (
                    f'{pheno_files_path_ct}{chrom}/{gene}_{celltype}_pheno_cov.tsv'
                )
                pheno_new_file = f'{condition_pheno_files_path}{celltype}/{chrom}/{gene}_{celltype}_conditional_pheno_cov.tsv'
                if not to_path(pheno_new_file).exists():
                    pheno_cond_job = get_batch().new_python_job(
                        name=f'gene make pheno cond file: {celltype},{gene}'
                    )
                    pheno_cond_job.storage(gene_condition_storage)
                    pheno_cond_job.memory(gene_condition_memory)
                    pheno_cond_job.call(
                        add_variant_to_pheno_file,
                        mt_path=chrom_mt_files_path,
                        gene=gene,
                        chrom=chrom,
                        celltype=celltype,
                        original_pheno_files_path=pheno_original_file,
                        new_pheno_files_path=to_path(pheno_new_file),
                        conditional_files_path=conditional_files_path,
                    )
                    manage_concurrency(pheno_cond_job)
                    logging.info(
                        f'make conditional pheno cov file job for {gene}, {celltype} scheduled'
                    )

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
