#!/usr/bin/env python3
# pylint: disable=import-error,line-too-long,too-many-arguments

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

analysis-runner \
    --config saige_assc.toml \
    --description "SAIGE-QTL association pipeline" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/output_files/" \
     python3 saige_assoc.py --celltypes CD4_Naive --chromosomes chr1 --vds-version vds1-0

"""

import click

import logging
import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, image_path, output_path


# Fit null model (step 1)
def build_fit_null_command(
    pheno_file: str,
    cov_col_list: str,
    sample_cov_col_list: str,
    sample_id_pheno: str,
    output_prefix: str,
    plink_path: str,
    pheno_col: str = 'y',
    trait_type: str = 'count',
    # this is a boolean but that's encoded differently between R and python
    skip_vre: str = 'FALSE',
    pheno_remove_zeros: str = 'FALSE',
    use_sparse_grm_null: str = 'FALSE',
    use_grm_null: str = 'FALSE',
    is_cov_offset: str = 'FALSE',
    is_cov_transform: str = 'TRUE',
    skip_model_fitting: str = 'FALSE',
    tol: float = 0.00001,
    is_overwrite_vre_file: str = 'TRUE',
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
    - trait type: count = Poisson, count_nb = Negative Binomial, quantitative = Normal
    - option to skip Variance Ratio estimation (discouraged)
    - option to add an offset to the fixed covariates (???)
    - option to transform (scale?) covariates?
    - option to skip model fitting (discouraged)
    - tolerance for convergence
    - overwrite variance ratio file (estimated here)

    Output:
    Rscript command (str) ready to run (bash)
    """
    pheno_file = get_batch().read_input(pheno_file)
    plink_prefix = get_batch().read_input_group(
        bim=f'{plink_path}.bim', bed=f'{plink_path}.bed', fam=f'{plink_path}.fam'
    )
    return f"""
        Rscript /usr/local/bin/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL={use_sparse_grm_null} \
        --useGRMtoFitNULL={use_grm_null} \
        --phenoFile={pheno_file} \
        --phenoCol={pheno_col} \
        --covarColList={cov_col_list} \
        --sampleCovarColList={sample_cov_col_list} \
        --sampleIDColinphenoFile={sample_id_pheno} \
        --traitType={trait_type} \
        --outputPrefix={output_prefix} \
        --skipVarianceRatioEstimation={skip_vre} \
        --isRemoveZerosinPheno={pheno_remove_zeros} \
        --isCovariateOffset={is_cov_offset} \
        --isCovariateTransform={is_cov_transform} \
        --skipModelFitting={skip_model_fitting} \
        --tol={tol} \
        --plinkFile={plink_prefix} \
        --IsOverwriteVarianceRatioFile={is_overwrite_vre_file}
    """


# Run single variant association (step 2)
def build_run_single_variant_test_command(
    output_path: str,
    vcf_file: str,
    chrom: str,
    cis_window_file: str,
    gmmat_model_path: str,
    variance_ratio_path: str,
    vcf_field: str = 'GT',
    min_maf: float = 0,
    min_mac: int = 5,
    loco_bool: str = 'FALSE',
    n_markers: int = 10000,
    spa_cutoff: int = 10000,
):
    """
    n.b. I have not edited this method docstring
            ... but I have definitely edited the method

    Build SAIGE command for running single variant test
    This will run a single variant test using a score test

    Input:
    - vcfFile: path to VCF file
    - vcfFileIndex: path to VCF index file (tbi)
    - saige output path: path to output saige file
    - chrom: chromosome to run this on
    - cis window: file with chrom | start | end to specify window
    - GMMAT model file: null model fit from previous step (.rda)
    - Variance Ratio file: as estimated from previous step (.txt)
    - group annotation: select only some annos from group file (e.g., lof)
    - group file: for each gene/set,
        one row specifying variants,
        one row specifying each variant's anno,
        one optional row with all weights
    - allele order: specifying whether alt-first or ref-first in genotype files
    - min MAF: minimum variant minor allele frequency to include
    - min MAC: minimum variant minor allele count to include
    - LOCO: leave one chromosome out (for what specifically?)
    - is no adjusted cov: covariate adjustment?
    - specify whether we're using a sparse GRM (vs full? vs no GRM at all?)

    Output:
    Rscript command (str) ready to run
    """

    if to_path(output_path).exists():
        return None, get_batch().read_input(output_path)

    vcf_group = get_batch().read_input_group(vcf=vcf_file, index=f'{vcf_file}.csi')
    second_job = get_batch().new_job(name="saige-qtl part 2")

    second_job.command(
        f"""
        Rscript /usr/local/bin/step2_tests_qtl.R \
        --vcfFile={vcf_group.vcf} \
        --vcfFileIndex={vcf_group.index} \
        --vcfField={vcf_field} \
        --SAIGEOutputFile={second_job.output} \
        --chrom={chrom} \
        --minMAF={min_maf} \
        --minMAC={min_mac} \
        --LOCO={loco_bool} \
        --GMMATmodelFile={gmmat_model_path} \
        --SPAcutoff={spa_cutoff} \
        --varianceRatioFile={variance_ratio_path} \
        --rangestoIncludeFile={cis_window_file} \
        --markers_per_chunk={n_markers}
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
        return None, get_batch().read_input(saige_gene_pval_output_file)

    saige_job = get_batch().new_job(name="saige-qtl part 3")
    saige_command_step3 = 'Rscript /usr/local/bin/step3_gene_pvalue_qtl.R'
    saige_command_step3 += f' --assocFile={saige_sv_output_file}'
    saige_command_step3 += f' --geneName={gene_name}'
    saige_command_step3 += f' --genePval_outputFile={saige_job.output}'
    saige_job.command(saige_command_step3)
    get_batch().write_output(saige_job.output, saige_gene_pval_output_file)
    return saige_job, saige_job.output


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
    cov_col_list: str,
    sample_cov_col_list: str,
    sample_id_pheno: str,
    plink_path: str,
    pheno_col: str,
):
    """
    Check if the output file already exists;
        if it does, read into a resource group
        if it doesn't, run the command and pass the resource group back
    Args:
        null_output_path ():
        pheno_file ():
        cov_col_list ():
        sample_cov_col_list ():
        sample_id_pheno ():
        plink_path ():
        pheno_col ():

    Returns:
        Tuple: (Job | None, ResourceGroup)

    """
    if to_path(null_output_path).exists():
        return None, get_batch().read_input_group(
            **{
                'rda': '{null_output_path}.rda',
                'varianceRatio': '{null_output_path}.varianceRatio.txt',
            }
        )
    gene_job = get_batch().new_job(name="saige-qtl")
    gene_job.image(image_path('saige-qtl'))
    apply_job_settings(gene_job, 'fit_null')

    # create output group for first command in gene job
    gene_job.declare_resource_group(
        output={
            'rda': '{root}.rda',
            'varianceRatio': '{root}.varianceRatio.txt',
        }
    )

    gene_job.command(
        build_fit_null_command(
            pheno_file=pheno_file,
            cov_col_list=cov_col_list,
            sample_cov_col_list=sample_cov_col_list,
            sample_id_pheno=sample_id_pheno,
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
    existing_cv_assoc_results = [
        file.name
        for file in to_path(gene_results_path).glob(f'{celltype}_*_cis_gene_pval')
    ]
    results_all_df = pd.concat(
        [
            pd.read_csv(to_path(pv_df), index_col=0)
            for pv_df in existing_cv_assoc_results
        ]
    )
    result_all_filename = to_path(out_path, 'analysis')
    logging.info(f'Write summary results to {result_all_filename}')
    with result_all_filename.open('w') as rf:
        results_all_df.to_csv(rf)


@click.command()
@click.option('--celltypes', help='add as one string, separated by comma')
@click.option('--chromosomes', help='add as one string, separated by comma')
@click.option('--vds-version', help=' e.g., 1-0 ')
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
@click.option('--sample-id', default='individual')
@click.option('--covs', default='sex,age,harmony_PC1,total_counts,sequencing_library')
@click.option('--sample-covs', default='sex,age')
@click.option('--cis-window-size', type=int, default=100000)
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=100,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
def main(
    celltypes: str,
    chromosomes: str,
    vds_version: str,
    # outputs from gene_expression processing
    pheno_cov_files_path: str,
    cis_window_files_path: str,
    # outputs from genotype processing
    genotype_files_prefix: str,
    # other
    sample_id: str,
    covs: str,
    sample_covs: str,
    cis_window_size: int,
    max_parallel_jobs: int = 100,
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
        if len(jobs) >= max_parallel_jobs:
            job.depends_on(jobs[-max_parallel_jobs])
        jobs.append(job)

    vre_plink_path = f'{genotype_files_prefix}/{vds_version}/vre_plink_2000_variants'

    for chromosome in chromosomes.split(','):

        # genotype vcf files are one per chromosome
        vcf_file_path = f'{genotype_files_prefix}/{vds_version}/{chromosome}_common_variants.vcf.bgz'
        # cis window files are split by gene but organised by chromosome also
        cis_window_files_path_chrom = f'{cis_window_files_path}/{chromosome}'

        for celltype in celltypes.split(','):

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

            # extract relevant gene-related files
            for gene in genes[0:2]:
                print(gene)
                pheno_cov_path = (
                    f'{pheno_cov_files_path_ct_chrom}/{gene}_{celltype}_pheno_cov.tsv'
                )
                cis_window_path = (
                    f'{cis_window_files_path_chrom}/{gene}_{cis_window_size}bp.tsv'
                )

                # check if these outputs already exist, if so don't make a new job
                null_job, null_output = run_fit_null_job(
                    output_path(f'output_files/{celltype}_{gene}'),
                    pheno_file=pheno_cov_path,
                    cov_col_list=covs,
                    sample_cov_col_list=sample_covs,
                    sample_id_pheno=sample_id,
                    plink_path=vre_plink_path,
                    pheno_col=gene,
                )
                if null_job:
                    manage_concurrency_for_job(null_job)

                # step 2 (cis eQTL single variant test)
                step2_job, step2_output = build_run_single_variant_test_command(
                    output_path=output_path(f'output_files/{celltype}_{gene}_cis'),
                    vcf_file=vcf_file_path,
                    chrom=chromosome,
                    cis_window_file=cis_window_path,
                    gmmat_model_path=null_output['rda'],
                    variance_ratio_path=null_output['varianceRatio'],
                )
                if step2_job:
                    manage_concurrency_for_job(step2_job)

                # step 3 (gene-level p-values)
                job3 = build_obtain_gene_level_pvals_command(
                    gene_name=gene,
                    saige_sv_output_file=step2_output,
                    saige_gene_pval_output_file=output_path(
                        f'output_files/{celltype}_{gene}_cis_gene_pval'
                    ),
                )

                if job3:
                    manage_concurrency_for_job(job3)

    # summarise results (per cell type)
    for celltype in celltypes.split(','):
        logging.info(f'start summarising results for {celltype}')
        summary_output_path = (
            f'output_files/summary_stats/{celltype}_all_cis_cv_results.tsv'
        )
        summarise_job = get_batch().new_python_job(
            f'Summarise CV results for {celltype}'
        )
        summarise_job.call(
            summarise_cv_results,
            gene_results_path=str('output_files/'),
            out_path=summary_output_path,
        )
        summarise_cv_results
    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
