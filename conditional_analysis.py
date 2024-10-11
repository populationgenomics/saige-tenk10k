#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script will run a conditional analysis

* Assumes a first complete run of SAIGE-QTL (common variants) has already run
    i.e., step1 (fit null), step2 (run association) and step3 (summarise)
* Runs step2 and step3 again, with the `condition` flag on
* requires
    * path to step1 results,
    * all step2 inputs as in the main pipeline (RV / CV)
    * one more file specifying gene-specific conditional strings

For multiple rounds of conditional analysis, this scripts will need to be re-run,
with an intermediate step that generates new conditional strings
This intermediate step is still TBD

To run (common variant test):

analysis-runner \
   --config saige_assoc_test.toml \
   --description "SAIGE-QTL CV conditional analysis" \
   --memory='32G' \
   --storage='50G' \
   --dataset "bioheart" \
   --access-level "full" \
   --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/output_files/240920" \
    python3 conditional_analysis.py --fit-null-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/output_files \
       --cis-window-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/input_files/cis_window_files \
       --genotype-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/input_files/genotypes/vds-tenk10k1-0 \
       --conditional-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/input_files/conditional_files \
       --common-or-rare='common'

To run (rare variant test):

analysis-runner \
   --config saige_assoc_test.toml \
   --description "SAIGE-QTL RV conditional analysis" \
   --memory='32G' \
   --storage='50G' \
   --dataset "bioheart" \
   --access-level "full" \
   --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/output_files/240920" \
    python3 conditional_analysis.py --fit-null-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/output_files \
       --group-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/input_files/group_files \
       --genotype-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/input_files/genotypes/vds-tenk10k1-0 \
       --conditional-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/240920/input_files/conditional_files \
       --common-or-rare='rare' --group-file-specs _dTSS_weights

"""

import click
import logging
import pandas as pd

from google.cloud import storage
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config, image_path, output_path
from cpg_utils.hail_batch import get_batch

# Run single variant or set-based association (step 2)
# with the `condition` flag
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


@click.option('--conditional-files-path', required=True)
@click.option('--fit-null-files-path', required=True)
@click.option('--genotype-files-path', required=True)
@click.option('--cis-window-or-group-files-path', required=True)
@click.option('--common-or-rare', default='common', help='type of analysis to perform')
@click.option('--group-file-specs', default='', help=' e.g. _dTSS_weights')
@click.option('--jobs-per-vm', type=int, default=25)
@click.command()
def conditional_analysis(
    # file path containing what snps to condition on for each gene
    conditional_files_path: str,
    # outputs from step 1 of saige
    fit_null_files_path: str,
    # cis window for common, group for rare
    cis_window_or_group_files_path: str,
    # outputs from get_genotype_vcf.py
    genotype_files_path: str,
    # whether to run a single-variant (for 'common' variants)
    # or set-based (for 'rare' variants) test
    common_or_rare: str,
    # if running rare test, weight types
    group_file_specs: str,
    # jobs to run per virtual machine
    jobs_per_vm: int,
):
    batch = get_batch('SAIGE-QTL conditional pipeline')
    jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(jobs) >= get_config()['saige']['max_parallel_jobs']:
            job.depends_on(jobs[-get_config()['saige']['max_parallel_jobs']])
        jobs.append(job)

    # pull principal args from config
    chromosomes: list[str] = get_config()['saige']['chromosomes']
    celltypes: list[str] = get_config()['saige']['celltypes']
    celltype_jobs: dict[str, list] = dict()
    cis_window_size = get_config()['saige']['cis_window_size']

    for chromosome in chromosomes:

        # genotype vcf files are one per chromosome
        vcf_file_path = f'{genotype_files_path}/{chromosome}_common_variants.vcf.bgz'

        # read in vcf file once per chromosome
        vcf_group = get_batch().read_input_group(
            vcf=vcf_file_path, index=f'{vcf_file_path}.csi'
        )
        step2_job = create_second_job(vcf_file_path)
        jobs_in_vm = 0

        # cis window and group files are split by gene but organised by chromosome also
        cis_window_or_group_files_path_chrom = (
            f'{cis_window_or_group_files_path}/{chromosome}'
        )

        step2_job = create_second_job(vcf_file_path)
        # jobs_in_vm = 0

        for celltype in celltypes:
            # extract gene list based on genes for which we have conditional files
            conditional_files_path_ct_file = (
                f'{conditional_files_path}/{celltype}_conditional_analysis.tsv'
            )
            conditional_df = pd.read_csv(
                output_path(conditional_files_path_ct_file), sep='\t'
            )

            genes = conditional_df['gene']
            logging.info(
                f'Preparing to run conditional analysis for these genes: {", ".join(genes)}'
            )

            # extract relevant gene-related files
            for gene in genes:

                gene_dependency = get_batch().new_job(f'Always run job for {gene}')
                gene_dependency.always_run()
                manage_concurrency_for_job(gene_dependency)

                conditional_string = conditional_df[conditional_df['gene'] == gene][
                    'variants_to_condition_on'
                ]
                # define gene-specific key
                suffix = conditional_files_path.split('/')[-1]
                test_key = f'{celltype}_{chromosome}_{gene}_{suffix}_{common_or_rare}'
                # define output
                test_output_path = output_path(test_key, 'analysis')

                # if the output exists, do nothing
                if to_path(f'{test_output_path}.set').exists():
                    continue

                # define gene-specific null output (rda + varianceRatio)
                null_path = (
                    f'{fit_null_files_path}/{celltype}/{chromosome}/{celltype}_{gene}'
                )
                # define gene-specific cis window or group file
                if common_or_rare == 'common':
                    cis_window_or_group_file = f'{cis_window_or_group_files_path_chrom}/{gene}_{cis_window_size}bp.tsv'
                elif common_or_rare == 'rare':
                    cis_window_or_group_file = f'{cis_window_or_group_files_path_chrom}/{gene}_{cis_window_size}bp{group_file_specs}.tsv'

                step2_output = build_run_conditional_analysis_command(
                    job=step2_job,
                    test_key=test_key,
                    test_output_path=test_output_path,
                    vcf_group=vcf_group,
                    chrom=(chromosome[3:]),
                    cis_window_or_group_file=cis_window_or_group_file,
                    gmmat_model_path=f'{null_path}.rda',
                    variance_ratio_path=f'{null_path}.varianceRatio.txt',
                    common_or_rare=common_or_rare,
                    conditional_string=conditional_string,
                )
                jobs_in_vm += 1
            if common_or_rare == 'common':
                # step 3 (gene-level p-values)
                saige_gene_pval_output_file = output_path(
                    f'{celltype}/{chromosome}/{celltype}_{gene}_cis_gene_pval'
                )
                if not to_path(saige_gene_pval_output_file).exists():
                    job3 = build_obtain_gene_level_pvals_command(
                        gene_name=gene,
                        saige_sv_output_file=step2_output,
                        saige_gene_pval_output_file=output_path(
                            f'{celltype}/{chromosome}/{celltype}_{gene}_cis_gene_pval'
                        ),
                    )
                    job3.depends_on(step2_job)
                    # add this job to the list of jobs for this cell type
                    celltype_jobs.setdefault(celltype, []).append(job3)

            if jobs_in_vm >= jobs_per_vm:
                step2_job = create_second_job(vcf_file_path)
                jobs_in_vm = 0

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    conditional_analysis()
