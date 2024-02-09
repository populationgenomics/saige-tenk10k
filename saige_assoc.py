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
     python3 saige_assoc.py --celltypes CD4_Naive --chromosomes chr1

"""

import click

# import glob

import hailtop.batch as hb

# import pandas as pd
from typing import List

# this may be part of production pipelines
# from cpg_workflows.utils import can_reuse
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
        --plinkFile={plink_path} \
        --IsOverwriteVarianceRatioFile={is_overwrite_vre_file}
    """


# Run single variant association (step 2)
def build_run_single_variant_test_command(
    vcf_file: str,
    vcf_file_index: str,
    vcf_field: str,
    saige_output_file: str,
    chrom: str,
    cis_window_file: str,
    gmmat_model_path: str,
    variance_ratio_path: str,
    min_maf: float = 0,
    min_mac: int = 5,
    loco_bool: str = 'FALSE',
    n_markers: int = 10000,
    spa_cutoff: int = 10000,
):
    """Build SAIGE command for running single variant test
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
    return f"""
        Rscript /usr/local/bin/step2_tests_qtl.R \
        --vcfFile={vcf_file} \
        --vcfFileIndex={vcf_file_index} \
        --vcfField={vcf_field} \
        --SAIGEOutputFile={saige_output_file} \
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


# Combine single variant associations at gene level (step 3)
def build_obtain_gene_level_pvals_command(
    gene_name: str,
    saige_sv_output_file: str,
    saige_gene_pval_output_file: str,
):
    """Build SAIGE command to obtain gene-level pvals
    Only for single-variant tests (Step3)
    combines single-variant p-values to obtain one gene
    level p-value

    Input:
    - ouput of previous step, association file (txt)
    - gene we need to aggregate results for (across SNPs)
    - path for output file
    """
    saige_command_step3 = 'Rscript /usr/local/bin/step3_gene_pvalue_qtl.R'
    saige_command_step3 += f' --assocFile={saige_sv_output_file}'
    saige_command_step3 += f' --geneName={gene_name}'
    saige_command_step3 += f' --genePval_outputFile={saige_gene_pval_output_file}'
    return saige_command_step3


def apply_job_settings(job, job_name: str):
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


# def run_fit_null_job(
#     b: hb.Batch,
#     null_output_path: str,
# ):
#     if can_reuse(null_output_path):
#         return None, b.read_input_group(
#             **{
#                 'rda': '{null_output_path}.rda',
#                 'varianceRatio': '{null_output_path}.varianceRatio.txt',
#             }
#         )
#     gene_job = b.new_job(name="saige-qtl")
#     gene_job.image(image_path('saige-qtl'))
#     apply_job_settings(gene_job, 'fit_null')

#     # create output group for first command in gene job
#     gene_job.declare_resource_group(
#         output={
#             'rda': '{root}.rda',
#             'varianceRatio': '{root}.varianceRatio.txt',
#         }
#     )

#     gene_job.command(
#         build_fit_null_command(
#             pheno_file=pheno_cov_filename,
#             cov_col_list=covs_list,
#             sample_cov_col_list=sample_covs_list,
#             sample_id_pheno=sample_id,
#             output_prefix=gene_job.output,
#             plink_path=plink_path,
#             pheno_col=gene_name,
#         )
#     )

#     # copy the output file to persistent storage
#     if null_output_path.startswith('gs://'):
#         b.write_output(gene_job.output, null_output_path)
#     return gene_job, gene_job.output


def association_pipeline(
    pheno_cov_filename: str,
    vcf_file_path: str,
    covs_list: str,
    sample_covs_list: str,
    sample_id: str,
    null_output_path: str,
    sv_output_path: str,
    gene_pvals_output_path: str,
    plink_path: str,
    gene_name: str,
    chrom: str,
    cis_window_file: str,
    vcf_field: str = 'GT',
):
    """
    Run association for one gene
    """

    batch = get_batch()

    gene_job = batch.new_job(name="saige-qtl")
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
            pheno_file=pheno_cov_filename,
            cov_col_list=covs_list,
            sample_cov_col_list=sample_covs_list,
            sample_id_pheno=sample_id,
            output_prefix=gene_job.output,
            plink_path=plink_path,
            pheno_col=gene_name,
        )
    )

    # copy the output file to persistent storage
    if null_output_path.startswith('gs://'):
        batch.write_output(gene_job.output, null_output_path)
    # gene_job, null_output_path = run_fit_null_job()
    # if gene_job:
    #     apply_job_settings(gene_job, 'fit_null')

    # step 2 (cis eQTL single variant test)
    vcf_group = batch.read_input_group(vcf=vcf_file_path, index=f'{vcf_file_path}.csi')
    gene_job.command(
        build_run_single_variant_test_command(
            vcf_file=vcf_group.vcf,
            vcf_file_index=vcf_group.index,
            vcf_field=vcf_field,
            saige_output_file=gene_job.output_single_variant,
            chrom=chrom,
            cis_window_file=cis_window_file,
            gmmat_model_path=gene_job.output['rda'],
            variance_ratio_path=gene_job.output['varianceRatio'],
        )
    )

    if sv_output_path.startswith('gs://'):
        batch.write_output(gene_job.output_single_variant, sv_output_path)

    # step3 (gene-level pvalues)
    gene_job.command(
        build_obtain_gene_level_pvals_command(
            gene_name=gene_name,
            saige_sv_output_file=gene_job.output_single_variant,
            saige_gene_pval_output_file=gene_job.get_gene_pvals,
        )
    )

    if gene_pvals_output_path.startswith('gs://'):
        batch.write_output(gene_job.get_gene_pvals, gene_pvals_output_path)


@click.command()
@click.option('--celltypes', help='add as one string, separated by comma')
@click.option('--chromosomes', help='add as one string, separated by comma')
@click.option(
    '--pheno-cov-files-path', default=output_path('input_files/pheno_cov_files/')
)
@click.option(
    '--cis-window-files-path', default=output_path('input_files/cis_window_files/')
)
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=100,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
@click.option('--cis-window-size', type=int, default=100000)
def main(
    celltypes: str,
    chromosomes: str,
    # outputs from gene_expression processing
    pheno_cov_files_path: str,
    cis_window_files_path: str,
    # outputs from genotype processing
    vcf_file_path: str,
    vre_plink_path: str,
    # other
    sample_id: str,
    covs: str,
    sample_covs: str,
    max_parallel_jobs: int = 100,
):
    """
    Run SAIGE-QTL pipeline for all cell types
    """

    batch = get_batch('SAIGE-QTL pipeline')
    jobs: List[hb.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(jobs) >= max_parallel_jobs:
            job.depends_on(jobs[-max_parallel_jobs])
        jobs.append(job)

    for chromosome in chromosomes.split(','):

        cis_window_files_path_chrom = f'{cis_window_files_path}/{chromosome}'

        for celltype in celltypes.split(','):

            # extract gene list based on genes for which we have pheno cov files
            pheno_cov_files_path_ct_chrom = (
                f'{pheno_cov_files_path}/{celltype}/{chromosome}'
            )
            files = to_path(pheno_cov_files_path_ct_chrom).glob(f'*_{celltype}.tsv')
            genes = [
                f.replace(pheno_cov_files_path_ct_chrom, '').replace(
                    f'_{celltype}.tsv', ''
                )
                for f in files
            ]

            # extract relevant gene-related files
            for gene in genes:
                pheno_cov_path = dataset_path(
                    f'{pheno_cov_files_path_ct_chrom}/{gene}_{celltype}.tsv'
                )
                cis_window_path = dataset_path(
                    f'{cis_window_files_path_chrom}/{gene}.tsv'
                )

                # todo - check if these outputs already exist, if so don't make a new job

                job = association_pipeline(
                    pheno_cov_filename=pheno_cov_path,
                    vcf_file_path=vcf_file_path,
                    covs_list=covs,
                    sample_covs_list=sample_covs,
                    sample_id=sample_id,
                    null_output_path=output_path(f'output_files/{celltype}_{gene}'),
                    sv_output_path=output_path(f'output_files/{celltype}_{gene}_cis'),
                    gene_pvals_output_path=output_path(
                        f'output_files/{celltype}_{gene}_cis_gene_pval'
                    ),
                    plink_path=vre_plink_path,
                    gene_name=gene,
                    chrom=chromosome,
                    cis_window_file=cis_window_path,
                )
                manage_concurrency_for_job(job)
    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
