#!/usr/bin/env python3

"""
Copy of https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/corr_matrix_maker.py
to rework for common variant SAIGE-QTL results

This script is used to create a correlation matrix between STR and SNP genotypes, which is required by fine-mapping methods.

Output: Correlation matrix, and list of variants (as appears in the matrix)

Workflow (for each cell type):
1) Extract genes where the eSTR is significant (default = FDR<0.05).
For each gene,
##### DON'T ACTUALLY DO THIS STEP2
2) Extract the STR and SNP coordinates where the association signal is p < 5e-4, to reduce computational burden of fine-mapper.
#########
3) Obtain the genotypes for each extracted STR and SNP in 2)
4) Calculate the correlation matrix between STR and SNP genotypes.

analysis-runner --dataset "bioheart" \
    --description "Calculate LD between STR and SNPs" \
    --access-level "test" \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr/fine_mapping/prep_files/v2" \
    corr_matrix_maker.py --snp-vcf-dir=gs://cpg-bioheart-test/str/associatr/tob_freeze_1/bgzip_tabix/v4 \
    --str-vcf-dir=gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1-chr-specific \
    --celltypes=gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive \
    --job-storage=10G \
    --max-parallel-jobs=50 \
    --str-fdr-dir=gs://cpg-bioheart-test-analysis/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat \
    --associatr-dir=gs://cpg-bioheart-test/str/associatr/snps_and_strs/rm_str_indels_dup_strs/tob_n1055_and_bioheart_n990/meta_results \
    --chromosomes=chr22


"""

import ast

import click
import numpy as np
import pandas as pd

import hailtop.batch as hb
from hailtop.batch import ResourceGroup

from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch


def ld_parser(
    snp_vcf_path: ResourceGroup,
    str_vcf_path: ResourceGroup,
    str_fdr: pd.DataFrame,
    celltype: str,
    pval_cutoff: float,
    associatr_dir: str,
) -> str:
    import pandas as pd
    from cyvcf2 import VCF

    if str_fdr.empty:
        print(f'No eSTRs for {celltype}')
        return  # type: ignore

    for index, row in str_fdr.iterrows():  # iterate over each gene
        gene = row['gene_name']
        chrom = ast.literal_eval(row['chr'])[0]
        # obtain raw associaTR results for this gene
        try:
            associatr_file = (
                f'{associatr_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv'
            )
            associatr = pd.read_csv(associatr_file, sep='\t')
        except FileNotFoundError:
            print(f'FileNotFound for this gene: {gene}')
            continue
        # apply p-value cutoff (mostly to reduce computational burden for fine-mapper later)
        associatr = associatr[associatr['pval_meta'] < pval_cutoff]
        if associatr.empty:
            print(f'No associatr results for this gene: {gene}')
            continue

        # create empty DF to store the relevant GTs (SNPs, STRs)
        snp_df = pd.DataFrame(columns=['individual'])
        str_df = pd.DataFrame(columns=['individual'])

        # cyVCF2 reads the SNP,STR VCF
        snp_vcf = VCF(snp_vcf_path['vcf'])
        snp_df['individual'] = snp_vcf.samples
        str_vcf = VCF(str_vcf_path['vcf'])
        snp_df['individual'] = snp_df['individual'].apply(
            lambda x: f"CPG{x}",
        )
        str_df['individual'] = str_vcf.samples
        str_df['individual'] = str_df['individual'].apply(
            lambda x: f"CPG{x}",
        )  # add CPG prefix to match SNP individual names

        for (
            index,
            row,
        ) in (
            associatr.iterrows()
        ):  # obtain GTs for each STR/SNP listed in the raw associatr file
            pos = row['pos']
            motif = row['motif']
            if '-' in row['motif']:  # SNP
                for variant in snp_vcf(f'{chrom}:{pos}-{pos}'):
                    if str(variant.INFO.get('RU')) == motif:
                        gt = variant.gt_types
                        gt[gt == 3] = 2  # convert HOM ALT 3 encoding into a 2

                        snp_df[f'{chrom}:{pos}_{motif}'] = gt
                        break
            else:  # STR
                chrom = row['chr']
                end = round(pos + row['ref_len'] * row['period'])
                for variant in str_vcf(f'{chrom}:{pos}-{pos}'):
                    if (str(variant.INFO.get('RU')) == motif) and (
                        int(variant.INFO.get('END')) == end
                    ):  # check if the motif and end coordinate matches
                        genotypes = variant.format('REPCN')
                        # Replace '.' with '-99/-99' to handle missing values
                        genotypes = np.where(genotypes == '.', '-99/-99', genotypes)

                        # Split each element by '/'
                        split_genotypes = [
                            genotype.split('/') for genotype in genotypes
                        ]

                        # Convert split_genotypes into a numpy array for easier manipulation
                        split_genotypes_array = np.array(split_genotypes)

                        # Convert the strings to integers and sum them row-wise
                        sums = np.sum(split_genotypes_array.astype(int), axis=1)
                        # set dummy -198 value to np.nan
                        sums = np.where(sums == -198, np.nan, sums)

                        str_df[f'{chrom}:{pos}_{motif}'] = sums
                        break
        # merge the two dataframes
        merged_df = str_df.merge(snp_df, on='individual')

        # calculate pairwise correlation of every variant
        merged_df = merged_df.drop(columns='individual')
        merged_df = merged_df.fillna(
            merged_df.mean(),
        )  # fill missing values with mean of the column (variant) to avoid NAs
        corr_matrix = merged_df.corr()

        # write to bucket
        corr_matrix.to_csv(
            output_path(
                f'correlation_matrix/{celltype}/{chrom}/{gene}_correlation_matrix.tsv',
                'analysis',
            ),
            sep='\t',
        )
        print(f"Wrote correlation matrix for {gene}")
    return  # type: ignore


@click.option(
    '--snp-vcf-dir',
    help='GCS file dir to SNP VCF files.',
    type=str,
)
@click.option(
    '--str-vcf-dir',
    help='GCS file dir to STR VCF files.',
    type=str,
)
@click.option(
    '--fdr-cutoff',
    help='FDR cutoff to use for eGene selection',
    type=float,
    default=0.05,
)
@click.option(
    '--pval-cutoff',
    help='P-value cutoff to use for associatr results (reduce finemapping computational burden)',
    default=5e-4,
)
@click.option('--celltypes', help='Cell types to use for coloc', type=str)
@click.option(
    '--str-fdr-dir',
    help='Path to STR FDR dir',
    default='gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat',
)
@click.option(
    '--associatr-dir',
    help='Path to STR-SNP associatr results',
    default='gs://cpg-bioheart-main-analysis/str/associatr/snps_and_strs/tob_n1055_and_bioheart_n990/meta_results',
)
@click.option(
    '--chromosomes',
    help='Chromosomes to use',
    default='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22',
)
@click.option('--job-cpu', default=1)
@click.option('--job-storage', default='20G')
@click.option('--max-parallel-jobs', default=22)
@click.command()
def main(
    snp_vcf_dir: str,
    str_vcf_dir: str,
    celltypes: str,
    fdr_cutoff: float,
    str_fdr_dir: str,
    job_cpu: int,
    job_storage: str,
    associatr_dir: str,
    pval_cutoff: float,
    chromosomes: str,
    max_parallel_jobs: int,
):
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    b = get_batch(name='Correlation matrix runner')
    for celltype in celltypes.split(','):
        # read in eSTR file
        str_fdr_file = f'{str_fdr_dir}/{celltype}_qval.tsv'
        str_fdr = pd.read_csv(str_fdr_file, sep='\t')
        str_fdr = str_fdr[
            str_fdr['qval'] < fdr_cutoff
        ]  # subset to eSTRs passing FDR 5% threshold by default
        for chrom in chromosomes.split(','):
            # filter eSTRs by chromosome
            str_fdr_chrom = str_fdr[str_fdr['chr'].str.contains("'" + chrom + "'")]
            # read in STR and SNP VCFs for this chromosome
            snp_vcf_path = f'{snp_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'
            str_vcf_path = f'{str_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'
            # run LD calculation for each chrom-celltype combination
            ld_job = b.new_python_job(
                f'LD calc for {celltype}:{chrom}',
            )
            ld_job.cpu(job_cpu)
            ld_job.storage(job_storage)
            snp_input = get_batch().read_input_group(
                **{'vcf': snp_vcf_path, 'tbi': snp_vcf_path + '.tbi'}
            )
            str_input = get_batch().read_input_group(
                **{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'}
            )

            ld_job.call(
                ld_parser,
                snp_input,
                str_input,
                str_fdr_chrom,
                celltype,
                pval_cutoff,
                associatr_dir,
            )
            manage_concurrency_for_job(ld_job)

    b.run(wait=False)


if __name__ == '__main__':
    main()
