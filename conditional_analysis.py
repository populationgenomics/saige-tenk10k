#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script will run a conditional analysis
"""

import click
import json
import pandas as pd

from google.cloud import storage
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config, image_path, output_path
from cpg_utils.hail_batch import get_batch

# Run single variant or set-based association (step 2)
# with condition
def build_run_conditional_analysis_command(
    job: hb.batch.job.Job,
    test_key: str,
    test_output_path: str,
    vcf_group: hb.ResourceGroup,
    chrom: str,
    cis_window_or_group_file: str,
    gmmat_model_path: str,
    variance_ratio_path: str,
    conditional_string: str,
    common_or_rare: str = 'common',
):
    """
    Build SAIGE command for running either a single variant test
    or a set-based test depending on the common_or_rare flag

    Input:
    - job: job to load this command into
    - test_key: unique key for this test
    - vcf_group: ResourceGroup with vcf and index
    - test output path: path to output saige file
    - chrom: chromosome to run this on
    - either / or
      - cis window: file with chrom | start | end to specify window
      - group: file with variants to test + weights + annotations
    - GMMAT model file: null model fit from previous step (.rda)
    - Variance Ratio file: as estimated from previous step (.txt)
    - SNPs to condition on (provided as a single comma-separated string)

    Output:
    Rscript command (str) ready to run
    """
    if common_or_rare == 'common':
        cis_window_file = get_batch().read_input(cis_window_or_group_file)
        variants_to_include_arg = f'--rangestoIncludeFile={cis_window_file}'
        args_from_config = ' '.join(
            [
                f'--{key}={value}'
                for key, value in get_config()['saige']['sv_test'].items()
            ]
        )
        # declare a uniquely named resource group for this single variant test
        job.declare_resource_group(**{test_key: {'output': f'{test_key}.output'}})
        output_arg = f'--SAIGEOutputFile={job[test_key].output}'
    elif common_or_rare == 'rare':
        group_file = get_batch().read_input(cis_window_or_group_file)
        variants_to_include_arg = f'--groupFile={group_file}'
        args_from_config = ' '.join(
            [
                f'--{key}={value}'
                for key, value in get_config()['saige']['set_test'].items()
            ]
        )
        # declare a uniquely named resource group for this set-based test
        rare_key_writeable = test_key.replace('/', '_')
        job.declare_resource_group(
            **{
                rare_key_writeable: {
                    'set': '{root}.set',
                    'singleAssoc.txt': '{root}.singleAssoc.txt',
                }
            }
        )
        output_arg = f'--SAIGEOutputFile={job[rare_key_writeable]}'

    job.command(
        f"""
        Rscript /usr/local/bin/step2_tests_qtl.R \
        --vcfFile={vcf_group.vcf} \
        --vcfFileIndex={vcf_group.index} \
        --chrom={chrom} \
        --GMMATmodelFile={gmmat_model_path} \
        --varianceRatioFile={variance_ratio_path} \
        --condition={conditional_string} \
        {variants_to_include_arg} \
        {output_arg} \
        {args_from_config}
    """
    )
    if common_or_rare == 'common':
        # write the output
        get_batch().write_output(job[test_key].output, test_output_path)
        return job[test_key].output
    elif common_or_rare == 'rare':
        get_batch().write_output(job[rare_key_writeable], test_output_path)


# Combine single variant associations at gene level (step 3)
def build_obtain_gene_level_pvals_command(
    gene_name: str,
    saige_sv_output_file: str,
    saige_gene_pval_output_file: str,
):
    """
    Build SAIGE command to obtain gene-level pvals
    Only for single-variant tests (Step3)
    combines single-variant p-values to obtain one gene
    level p-value

    Input:
    - output of previous step, association file (txt)
    - gene we need to aggregate results for (across SNPs)
    - path for output file
    """
    saige_job = get_batch().new_job(name="saige-qtl part 3")
    saige_command_step3 = f"""
        Rscript /usr/local/bin/step3_gene_pvalue_qtl.R \
        --assocFile={saige_sv_output_file} \
        --geneName={gene_name} \
        --genePval_outputFile={saige_job.output}
    """
    saige_job.image(image_path('saige-qtl'))
    saige_job.command(saige_command_step3)
    get_batch().write_output(saige_job.output, saige_gene_pval_output_file)
    return saige_job


def apply_job_settings(job: hb.batch.job.Job, job_name: str):
    """
    Apply settings to a job based on the name

    Args:
        job (hb.batch.job): job to apply settings to
        job_name (str): name used to find settings
    """
    # if this job has a section in config
    if job_settings := get_config()['saige']['job_specs'].get(job_name):
        # if memory setting in config - apply it
        if memory := job_settings.get('memory'):
            job.memory(memory)
        # if storage setting in config - apply it
        if storage := job_settings.get('storage'):
            job.storage(storage)
        # if cpu setting in config - apply it
        if cpu := job_settings.get('cpu'):
            job.cpu(cpu)


def create_second_job(vcf_path: str) -> hb.batch.job.Job:
    """
    Create a second job to run the single variant test
    """
    # get the size of the vcf file
    storage_client = storage.Client()
    bucket, filepath = vcf_path.removeprefix('gs://').split('/', 1)
    blob = storage_client.bucket(bucket).blob(filepath)
    blob.reload()  # refresh the blob to get the metadata
    size = blob.size // (1024**3)  # bytes to GB

    second_job = get_batch().new_job(name="saige-qtl part 2")
    apply_job_settings(second_job, 'sv_test')

    # VCF size, plus a 5GB buffer
    second_job.storage(f'{size + 5 }Gi')
    second_job.image(image_path('saige-qtl'))
    return second_job

def conditional_analysis_until_significant_round(
        gene_dict,
        significance_fdr_threshold=0.05,
    ):
    """
    Run conditional analysis for a gene + celltype combo
    continue to add SNPs to conditional analysis
    while there is a significant signal
    """
    # run step2
    # run step3
    # run fdr (how?)
    # check fdr<threshold
    # if any snp is significant, take top snp add to condition
    # repeat

@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--round1-results-path', required=True)
@click.option('--fit-null-files-path', required=True)
@click.option('--genotype-files-prefix', required=True)
@click.option('--cis-window-or-group-files-path', required=True)
@click.option('--qv-significance-threshold', default=0.05)
@click.option('--common-or-rare', default='common', help='type of analysis to perform')
@click.option('--cis-window-size', default=100000)
@click.option('--group-file-specs', default='_dTSS')
@click.command()
def conditional_analysis(
    celltypes: str,
    chromosomes: str,
    # results from running saige-qtl main pipeline
    round1_results_path: str,
    # outputs from step 1 of saige
    fit_null_files_path: str,
    # cis window for common, group for rare
    cis_window_or_group_files_path: str,
    # outputs from get_genotype_vcf.py
    genotype_files_prefix: str,
    # significance threshold
    qv_significance_threshold: float,
    common_or_rare: str,
    cis_window_size: int,
    group_file_specs: str,
):
    batch = get_batch('SAIGE-QTL conditional pipeline')

    # for every cell type, have a dictionary that contains all genes with a significant eQTL
    # and info on the top SNP
    # this will get updated by progressive rounds of conditional analysis
    # so every gene will get (or not) one more SNP added
    # should this be a table instead with columns round 1, round 2 etc?
    for celltype in celltypes.split(','):
        # open the summary results from round 1 of running SAIGE-QTL (common variant analysis)
        celltype_results_path = f'{round1_results_path}/summary_stats/{celltype}_common_top_snp_cis_raw_pvalues.tsv'
        results_df = pd.read_csv(celltype_results_path)
        # extract significant results
        results_df_sign = results_df[results_df['qv'] < qv_significance_threshold]
        genes = results_df_sign[['gene']]
        top_snps = results_df_sign[['top_snp']]
        # as more rounds of conditional analysis are performed, more snps will be added
        significant_snps_gene_dict: dict = {genes: top_snps}
        # temporarily write this out?
        significant_genes_dict_path = output_path(
            f'conditional_analysis/{celltype}/round1_significant_genes.json'
        )
        with to_path(significant_genes_dict_path).open('wt') as out_handle:
            json.dump(significant_genes_dict_path, fp=out_handle, indent=4)

    # next, for every chromosome (because the VCF are saved chromosome-wise) and cell type
    # get the genes with at least a significant SNP from the dicts above
    # look up the intermediate files from step 1 of SAIGE-QTL (fit null)
    # run step 2 with the top SNP as condition

    # Q1: how do I do this iteratively for consequent conditional rounds?
    # Q2: how do I separate from the ability to just run a conditional analysis in a more ad hoc way?
    # for example, a) conditioning on common eQTLs when running rare variant test
    # b) conditioning on top eQTL from cell type A when running cell type B, to assess cell type specificity
    for chromosome in chromosomes:

        # genotype vcf files are one per chromosome
        vcf_file_path = f'{genotype_files_prefix}/{chromosome}_common_variants.vcf.bgz'

        # read in vcf file once per chromosome
        vcf_group = get_batch().read_input_group(
            vcf=vcf_file_path, index=f'{vcf_file_path}.csi'
        )

        # cis window and group files are split by gene but organised by chromosome also
        cis_window_or_group_files_path_chrom = (
            f'{cis_window_or_group_files_path}/{chromosome}'
        )

        step2_job = create_second_job(vcf_file_path)
        # jobs_in_vm = 0

        for celltype in celltypes:
            # extract significant genes from json
            for gene in genes:
                # define gene-specific key
                test_key = f'{celltype}_{chromosome}_{celltype}_{gene}_conditional_round{conditional_round_number}'
                # define output
                test_output_path = 'xxx'
                # define gene-specific null output (rda + varianceRatio)
                null_path = (
                    f'{fit_null_files_path}/{celltype}/{chromosome}/{celltype}_{gene}'
                )
                # define gene-specific cis window or group file
                if common_or_rare == 'common':
                    cis_window_or_group_file = f'{cis_window_or_group_files_path_chrom}/{gene}_{cis_window_size}bp.tsv'
                elif common_or_rare == 'rare':
                    cis_window_or_group_file = f'{cis_window_or_group_files_path_chrom}/{gene}_{cis_window_size}bp{group_file_specs}.tsv'
                # build command
                # add to job

                build_run_conditional_analysis_command(
                    job=step2_job,
                    test_key=test_key,
                    rare_output_path=test_output_path,
                    vcf_group=vcf_group,
                    chrom=(chromosome[3:]),
                    cis_window_or_group_file=cis_window_or_group_file,
                    gmmat_model_path=f'{null_path}.rda',
                    variance_ratio_path=f'{null_path}.varianceRatio.txt',
                )


if __name__ == '__main__':
    conditional_analysis()
