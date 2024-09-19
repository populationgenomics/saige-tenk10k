#!/usr/bin/env python3
# pylint: disable=import-error,line-too-long

"""
Hail Batch workflow to perform association tests using SAIGE-QTL

- given all input files already generated
    - pheno + cov file (from get_anndata.py)
    - cis window file (from get_anndata.py)
    - genotype file (from get_genotype_vcf.py)
    - VRE genotypes (from get_genotype_vcf.py)
- builds saige commands (just add to str)
- run SAIGE-QTL (execute Rscript from command line)
- aggregate & summarise results (not yet)


To run:

In test:

analysis-runner \
    --config saige_assoc_test.toml \
    --description "SAIGE-QTL association pipeline" \
    --memory='32G' \
    --storage='50G' \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/v1" \
     python3 saige_assoc.py \
     --pheno-cov-files-path=gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/input_files/pheno_cov_files \
        --cis-window-files-path=gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/input_files/cis_window_files \
        --genotype-files-prefix=gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/input_files/genotypes/v3/vds-tenk10k1-0 \
        --vre-files-prefix=gs://cpg-bioheart-test/saige-qtl/bioheart_n990_and_tob_n1055/input_files/genotypes/v3/vds-tenk10k1-0


In main:

analysis-runner \
    --config saige_assoc_test.toml \
    --description "SAIGE-QTL association pipeline" \
    --memory='32G' \
    --storage='50G' \
    --dataset "bioheart" \
    --access-level "full" \
    --output-dir "saige-qtl/bioheart_n990/v1" \
     python3 saige_assoc.py \
     --pheno-cov-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990/input_files/v1/pheno_cov_files \
        --cis-window-files-path=gs://cpg-bioheart-main/saige-qtl/bioheart_n990/input_files/v1/cis_window_files \
        --genotype-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990/input_files/genotypes/vds-bioheart1-0 \
        --vre-files-prefix=gs://cpg-bioheart-main/saige-qtl/bioheart_n990/input_files/genotypes/vds-bioheart1-0


"""

import click

import logging

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, image_path, output_path


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


# Run single variant association (step 2)
def build_run_single_variant_test_command(
    output_path: str,
    vcf_file: str,
    chrom: str,
    cis_window_file: str,
    gmmat_model_path: str,
    variance_ratio_path: str,
):
    """
    Build SAIGE command for running single variant test
    This will run a single variant test using a score test

    Input:
    - vcfile: path to VCF file
    - vcfFileIndex: path to VCF index file (csi)
    - saige output path: path to output saige file
    - chrom: chromosome to run this on
    - cis window: file with chrom | start | end to specify window
    - GMMAT model file: null model fit from previous step (.rda)
    - Variance Ratio file: as estimated from previous step (.txt)

    Output:
    Rscript command (str) ready to run
    """

    if to_path(output_path).exists():
        return None, get_batch().read_input(output_path)

    vcf_group = get_batch().read_input_group(vcf=vcf_file, index=f'{vcf_file}.csi')
    cis_window_file = get_batch().read_input(cis_window_file)

    second_job = get_batch().new_job(name="saige-qtl part 2")
    apply_job_settings(second_job, 'sv_test')
    second_job.image(image_path('saige-qtl'))

    args_from_config = ' '.join(
        [f'--{key}={value}' for key, value in get_config()['saige']['sv_test'].items()]
    )

    second_job.command(
        f"""
        Rscript /usr/local/bin/step2_tests_qtl.R \
        --vcfFile={vcf_group.vcf} \
        --vcfFileIndex={vcf_group.index} \
        --SAIGEOutputFile={second_job.output} \
        --chrom={chrom} \
        --GMMATmodelFile={gmmat_model_path} \
        --varianceRatioFile={variance_ratio_path} \
        --rangestoIncludeFile={cis_window_file} \
        {args_from_config}
    """
    )

    # write the output
    get_batch().write_output(second_job.output, output_path)

    return second_job, second_job.output


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
    if to_path(saige_gene_pval_output_file).exists():
        return None

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


def summarise_cv_results(
    celltype: str,
    gene_results_path: str,
    out_path: str,
):
    """
    Summarise gene-specific results
    """
    import logging
    import pandas as pd
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    existing_cv_assoc_results = [
        str(file)
        for file in to_path(gene_results_path).glob(f'*/{celltype}_*_cis_gene_pval')
    ]
    results_all_df = pd.concat(
        [
            pd.read_csv(to_path(pv_df), index_col=0)
            for pv_df in existing_cv_assoc_results
        ]
    )
    result_all_filename = to_path(output_path(out_path, category='analysis'))
    logging.info(f'Write summary results to {result_all_filename}')
    with result_all_filename.open('w') as rf:
        results_all_df.to_csv(rf)


@click.option(
    '--pheno-cov-files-path',
    default=dataset_path('saige-qtl/input_files/pheno_cov_files'),
)
@click.option(
    '--cis-window-files-path',
    default=dataset_path('saige-qtl/input_files/cis_window_files'),
)
@click.option(
    '--genotype-files-prefix', default=dataset_path('saige-qtl/input_files/genotypes')
)
@click.option(
    '--vre-files-prefix', default=dataset_path('saige-qtl/input_files/genotypes')
)
@click.option('--test-str', is_flag=True, help='Test with STR VCFs')
@click.command()
def main(
    # outputs from gene_expression processing
    pheno_cov_files_path: str,
    cis_window_files_path: str,
    # outputs from genotype processing
    genotype_files_prefix: str,
    vre_files_prefix: str,
    test_str: bool = False,
):
    """
    Run SAIGE-QTL pipeline for all cell types
    """

    batch = get_batch('SAIGE-QTL pipeline')
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

    vre_plink_path = f'{vre_files_prefix}/vre_plink_2000_variants'

    for chromosome in chromosomes:
        # genotype vcf files are one per chromosome
        if test_str:
            vcf_file_path = (
                f'{genotype_files_prefix}/hail_filtered_{chromosome}.vcf.bgz'
            )
        else:
            vcf_file_path = (
                f'{genotype_files_prefix}/{chromosome}_common_variants.vcf.bgz'
            )
        # cis window files are split by gene but organised by chromosome also
        cis_window_files_path_chrom = f'{cis_window_files_path}/{chromosome}'

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
                cis_window_path = (
                    f'{cis_window_files_path_chrom}/{gene}_{cis_window_size}bp.tsv'
                )

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

                # step 2 (cis eQTL single variant test)
                step2_job, step2_output = build_run_single_variant_test_command(
                    output_path=output_path(
                        f'{celltype}/{chromosome}/{celltype}_{gene}_cis', 'analysis'
                    ),
                    vcf_file=vcf_file_path,
                    chrom=(chromosome[3:]),
                    cis_window_file=cis_window_path,
                    gmmat_model_path=null_output['rda'],
                    variance_ratio_path=null_output['varianceRatio.txt'],
                )
                if step2_job is not None:
                    step2_job.depends_on(gene_dependency)
                    gene_dependency = step2_job

                # step 3 (gene-level p-values)
                job3 = build_obtain_gene_level_pvals_command(
                    gene_name=gene,
                    saige_sv_output_file=step2_output,
                    saige_gene_pval_output_file=output_path(
                        f'{celltype}/{chromosome}/{celltype}_{gene}_cis_gene_pval'
                    ),
                )

                if job3 is not None:
                    job3.depends_on(gene_dependency)
                    gene_dependency = job3

                # add this job to the list of jobs for this cell type
                if gene_dependency is not None:
                    celltype_jobs.setdefault(celltype, []).append(gene_dependency)

    # summarise results (per cell type)
    for celltype in celltypes:
        logging.info(f'start summarising results for {celltype}')
        summary_output_path = (
            f'summary_stats/{celltype}_all_cis_cv_gene_level_results.tsv'
        )

        summarise_job = get_batch().new_python_job(
            f'Summarise CV results for {celltype}'
        )
        summarise_job.depends_on(*celltype_jobs[celltype])
        summarise_job.call(
            summarise_cv_results,
            celltype=celltype,
            gene_results_path=output_path(celltype),
            out_path=summary_output_path,
        )
    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
