#!/usr/bin/env python3
# pylint: disable=import-error,line-too-long

"""
Hail Batch workflow to perform set-based tests using SAIGE-QTL

Similar to saige_assoc.py but for set-based tests
typically used for rare variants

- given all input files already generated
    - pheno + cov file (from get_anndata.py)
    - genotype file (from get_genotype_vcf.py)
    - VRE genotypes (from get_genotype_vcf.py)
    - group file (from make_group_file.py)
- builds saige commands (just add to str)
- run set-based tests using SAIGE-QTL (execute Rscript from command line)
- aggregate & summarise results


To run:

analysis-runner \
   --config saige_assoc_test.toml \
   --description "SAIGE-QTL RV association pipeline" \
   --memory='32G' \
   --storage='50G' \
   --dataset "bioheart" \
   --access-level "full" \
   --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/output_files/240920" \
    python3 saige_assoc_set_test.py \
    --pheno-cov-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/pheno_cov_files \
       --group-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/group_files \
       --genotype-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/genotypes/vds-tenk10k1-0 \
       --vre-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920/genotypes/vds-tenk10k1-0 \
       --group-file-specs _dTSS_weights

"""

import click
import json
import logging

from datetime import datetime
from os import getenv

from google.cloud import storage
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config, image_path, output_path, try_get_ar_guid
from cpg_utils.hail_batch import dataset_path, get_batch


# Fit null model (step 1)
def build_fit_null_command(
    pheno_file: str, output_prefix: str, plink_path: str, pheno_col: str
):
    """Build SAIGE command for fitting null model
    This will fit a Poisson / NB mixed model under the null hypothesis

    Input:
    - Phenotype / covariate file - rows: samples, cols: pheno (y), cov1, cov2 etc
    - Comma separated str of column names from previous file to be used as covariates
    - Same as above but for sample specific covariates
    - Column name specifying sample / individual (e.g., IND_ID, or CPG_ID etc)
    - output prefix: where to save the fitted model (.rda)
    - Plink path: path to plink file (subset of ~2,000 markers for VRE)
    - pheno col: name of column specifying pheno (default: "y")
    - overwrite variance ratio file (estimated here)

    Output:
    Rscript command (str) ready to run (bash)
    """
    pheno_file = get_batch().read_input(pheno_file)
    plink_prefix = get_batch().read_input_group(
        bim=f'{plink_path}.bim', bed=f'{plink_path}.bed', fam=f'{plink_path}.fam'
    )

    # pull all values from the config file's saige.build_fit_null section
    args_from_config = ' '.join(
        [
            f'--{key}={value}'
            for key, value in get_config()['saige']['build_fit_null'].items()
        ]
    )

    return f"""
        Rscript /usr/local/bin/step1_fitNULLGLMM_qtl.R \
        --phenoFile={pheno_file} \
        --plinkFile={plink_prefix} \
        --outputPrefix={output_prefix} \
        --phenoCol={pheno_col} \
        {args_from_config }
    """


# Run set-based association (step 2b)
def build_run_set_based_test_command(
    job: hb.batch.job.Job,
    rare_key: str,
    rare_output_path: str,
    vcf_group: hb.ResourceGroup,
    chrom: str,
    group_file: str,
    gmmat_model_path: str,
    variance_ratio_path: str,
):
    """
    Build SAIGE command for running set-based test
    This will run a single variant test using a score test

    Input:
    - job: job to run this command in
    - rare_key: unique key for this rare variant test (used to name output)
    - vcf group: VCF & index file ResourceGroup
    - set test output path: path to output saige file
    - chrom: chromosome to run this on
    - group file: file with variants to test, and weights
    - GMMAT model file: null model fit from previous step (.rda)
    - Variance Ratio file: as estimated from previous step (.txt)

    Output:
    Rscript command (str) ready to run
    """

    group_file = get_batch().read_input(group_file)

    args_from_config = ' '.join(
        [f'--{key}={value}' for key, value in get_config()['saige']['set_test'].items()]
    )

    # declare a uniquely named resource group for this set-based test
    rare_key_writeable = rare_key.replace('/', '_')
    job.declare_resource_group(
        **{
            rare_key_writeable: {
                'set': '{root}.set',
                'singleAssoc.txt': '{root}.singleAssoc.txt',
            }
        }
    )

    job.command(
        f"""
        Rscript /usr/local/bin/step2_tests_qtl.R \
        --vcfFile={vcf_group.vcf} \
        --vcfFileIndex={vcf_group.index} \
        --SAIGEOutputFile={job[rare_key_writeable]} \
        --chrom={chrom} \
        --GMMATmodelFile={gmmat_model_path} \
        --varianceRatioFile={variance_ratio_path} \
        --groupFile={group_file} \
        {args_from_config}
    """
    )

    # write the output
    get_batch().write_output(job[rare_key_writeable], rare_output_path)


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


def run_fit_null_job(
    null_output_path: str,
    pheno_file: str,
    plink_path: str,
    pheno_col: str = 'y',
):
    """
    Check if the output file already exists;
        if it does, read into a resource group
        if it doesn't, run the command and pass the resource group back
    Args:
        null_output_path ():
        pheno_file ():
        plink_path ():
        pheno_col (str): defaults to "y", or supplied with a gene name

    Returns:
        Tuple: (Job | None, ResourceGroup)

    """
    if to_path(f'{null_output_path}.rda').exists():
        return None, get_batch().read_input_group(
            **{
                'rda': f'{null_output_path}.rda',
                'varianceRatio.txt': f'{null_output_path}.varianceRatio.txt',
            }
        )

    gene_job = get_batch().new_job(name="saige-qtl part 1")
    gene_job.image(image_path('saige-qtl'))
    apply_job_settings(gene_job, 'fit_null')

    # create output group for first command in gene job
    gene_job.declare_resource_group(
        output={
            'rda': '{root}.rda',
            'varianceRatio.txt': '{root}.varianceRatio.txt',
        }
    )

    gene_job.command(
        build_fit_null_command(
            pheno_file=pheno_file,
            output_prefix=gene_job.output,
            plink_path=plink_path,
            pheno_col=pheno_col,
        )
    )

    # copy the output file to persistent storage
    if null_output_path:
        get_batch().write_output(gene_job.output, null_output_path)

    return gene_job, gene_job.output


def summarise_rv_results(
    celltype: str,
    gene_results_path: str,
    summary_output_path: str,
):
    """
    Summarise gene-specific results
    """
    import logging
    import pandas as pd
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    existing_rv_assoc_results = [
        str(file)
        for file in to_path(gene_results_path).glob(f'*/{celltype}_*_cis_rare.set')
    ]
    results_all_df = pd.concat(
        [
            pd.read_csv(to_path(pv_df), index_col=0)
            for pv_df in existing_rv_assoc_results
        ]
    )
    result_all_filename = to_path(output_path(summary_output_path, category='analysis'))
    logging.info(f'Write summary results to {result_all_filename}')
    with result_all_filename.open('w') as rf:
        results_all_df.to_csv(rf)


def create_a_2b_job() -> hb.batch.job.Job:
    """
    Create a job that will run a set-based test
    """
    second_job = get_batch().new_job(name="saige-qtl part 2b")
    apply_job_settings(second_job, 'set_test')
    second_job.image(image_path('saige-qtl'))
    return second_job


@click.option(
    '--pheno-cov-files-path',
    default=dataset_path('saige-qtl/input_files/pheno_cov_files'),
)
@click.option(
    '--group-files-path',
    default=dataset_path('saige-qtl/input_files/group_files'),
)
@click.option(
    '--genotype-files-prefix', default=dataset_path('saige-qtl/input_files/genotypes')
)
@click.option(
    '--vre-files-prefix', default=dataset_path('saige-qtl/input_files/genotypes')
)
@click.option('--writeout-file-prefix', default=dataset_path('saige-qtl'))
@click.option('--ngenes-to-test', default='all')
@click.option('--group-file-specs', default='')
@click.option('--jobs-per-vm', default=10, type=int)
@click.command()
def main(
    # output from get_anndata.py
    pheno_cov_files_path: str,
    # output from make_group_file.py
    group_files_path: str,
    # outputs from get_genotype_vcf.py
    genotype_files_prefix: str,
    vre_files_prefix: str,
    # write out inputs and flags used for this run
    writeout_file_prefix: str,
    ngenes_to_test: str,
    group_file_specs: str,
    jobs_per_vm: int,
):
    """
    Run SAIGE-QTL RV pipeline for all cell types
    """

    batch = get_batch('SAIGE-QTL RV pipeline')
    jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(jobs) >= get_config()['saige']['max_parallel_jobs']:
            job.depends_on(jobs[-get_config()['saige']['max_parallel_jobs']])
        jobs.append(job)

    # define writeout file by type of pipeline and date and time
    date_and_time = datetime.today().strftime('%Y-%m-%d_%H:%M:%S')
    writeout_file = (
        f'{writeout_file_prefix}/saige_qtl_rare_variant_pipeline_{date_and_time}.csv'
    )

    # pull principal args from config
    chromosomes: list[str] = get_config()['saige']['chromosomes']
    celltypes: list[str] = get_config()['saige']['celltypes']
    celltype_jobs: dict[str, list] = dict()

    vre_plink_path = f'{vre_files_prefix}/vre_plink_2000_variants'

    # populate all the important params into a file for long-term reference
    writeout_dict: dict = {
        'ar_guid': try_get_ar_guid() or 'UNKNOWN',
        'results_output_path': output_path(''),
        'vre_plink_files_prefix_used': vre_plink_path,
        'pheno_cov_files_path_used': pheno_cov_files_path,
        'group_files_path_used': group_files_path,
        'saige_fit_null_params': get_config()['saige'],
        'saige_set_test_params': get_config()['saige.set_test'],
        'saige_fit_null_job_specs': get_config()['saige.job_specs.fit_null'],
        'saige_set_test_job_specs': get_config()['saige.job_specs.set_test'],
        'runtime_config': getenv('CPG_CONFIG_PATH'),
    }

    for chromosome in chromosomes:
        # create one job, and stack multiple jobs on it
        step2_job = create_a_2b_job()

        # to start with, we have no jobs in the image
        jobs_in_vm = 0
        # genotype vcf files are one per chromosome, so read in at the top
        vcf_file_path = f'{genotype_files_prefix}/{chromosome}_rare_variants.vcf.bgz'
        vcf_group = get_batch().read_input_group(
            vcf=vcf_file_path, index=f'{vcf_file_path}.csi'
        )
        writeout_dict[f'{chromosome}_vcf_file_used'] = vcf_file_path

        # group files are split by gene but organised by chromosome also
        group_files_path_chrom = f'{group_files_path}/{chromosome}'

        cis_window_size = get_config()['saige']['cis_window_size']

        for celltype in celltypes:
            # extract gene list based on genes for which we have pheno cov files
            pheno_cov_files_path_ct_chrom = (
                f'{pheno_cov_files_path}/{celltype}/{chromosome}'
            )
            logging.info(f'globbing {pheno_cov_files_path_ct_chrom}')

            # do a glob, then pull out all file names as Strings
            files = [
                file.name
                for file in to_path(pheno_cov_files_path_ct_chrom).glob(
                    f'*_{celltype}_pheno_cov.tsv'
                )
            ]
            # if specified, only test ngenes genes
            if ngenes_to_test != 'all':
                files = files[0 : int(ngenes_to_test)]
            logging.info(f'I found these files: {", ".join(files)}')

            genes = [f.replace(f'_{celltype}_pheno_cov.tsv', '') for f in files]
            logging.info(f'I found these genes: {", ".join(genes)}')

            drop_genes: list[str] = get_config()['saige']['drop_genes']

            genes = [x for x in genes if x not in drop_genes]

            # extract relevant gene-related files
            for gene in genes:
                pheno_cov_path = (
                    f'{pheno_cov_files_path_ct_chrom}/{gene}_{celltype}_pheno_cov.tsv'
                )
                group_path = f'{group_files_path_chrom}/{gene}_{cis_window_size}bp{group_file_specs}.tsv'

                gene_dependency = get_batch().new_job(f' Always run job for {gene}')
                gene_dependency.always_run()
                manage_concurrency_for_job(gene_dependency)

                # check if these outputs already exist, if so don't make a new job
                null_job, null_output = run_fit_null_job(
                    output_path(f'{celltype}/{chromosome}/{celltype}_{gene}'),
                    pheno_file=pheno_cov_path,
                    plink_path=vre_plink_path,
                    pheno_col=gene,
                )
                if null_job is not None:
                    null_job.depends_on(gene_dependency)
                    gene_dependency = null_job

                # step 2 (cis eQTL set-based test)
                # unique key for this set-based test
                rare_key = f'{celltype}/{chromosome}/{celltype}_{gene}_cis_rare'
                # unique output path for this set-based test
                rare_output_path = output_path(rare_key, 'analysis')

                # if the output exists, do nothing
                if to_path(f'{rare_output_path}.set').exists():
                    continue

                # instruct an additional command to run inside this VM
                build_run_set_based_test_command(
                    job=step2_job,
                    rare_key=rare_key,
                    rare_output_path=rare_output_path,
                    vcf_group=vcf_group,
                    chrom=(chromosome[3:]),
                    group_file=group_path,
                    gmmat_model_path=null_output['rda'],
                    variance_ratio_path=null_output['varianceRatio.txt'],
                )
                jobs_in_vm += 1

                # if we ran this job, an additional dependency
                step2_job.depends_on(gene_dependency)

                # add this job to the list of jobs for this cell type
                celltype_jobs.setdefault(celltype, []).append(step2_job)

                # check if we need a new VM, i.e. we've hit the jobs-per-VM limit
                if jobs_in_vm >= jobs_per_vm:
                    step2_job = create_a_2b_job()
                    jobs_in_vm = 0

    # summarise results (per cell type)
    for celltype in celltypes:
        logging.info(f'start summarising results for {celltype}')
        summary_output_path = (
            f'summary_stats/{celltype}_all_cis_rv_set_test_results.tsv'
        )

        summarise_job = get_batch().new_python_job(
            f'Summarise RV results for {celltype}'
        )
        summarise_job.depends_on(*celltype_jobs[celltype])
        summarise_job.call(
            summarise_rv_results,
            celltype=celltype,
            gene_results_path=output_path(celltype, category='analysis'),
            summary_output_path=summary_output_path,
        )

    # write the file containing all important input parameters
    with to_path(writeout_file).open('wt') as out_handle:
        json.dump(writeout_dict, fp=out_handle, indent=4)

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
