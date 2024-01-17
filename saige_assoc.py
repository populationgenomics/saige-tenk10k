#!/usr/bin/env python3
# pylint: disable=import-error,line-too-long,too-many-arguments

"""
Hail Batch workflow to perform association tests using SAIGE-QTL

- given all input files already generated
    - pheno_cov filename, genotypes with annotations
- builds saige commands (just add to str)
- run saige commands (execute Rscript from command line)
- aggregate & summarise results (not yet)


To run:

analysis-runner \
    --description "SAIGE-QTL association pipeline" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/repo-example-inputs/output/" \
     python3 saige_assoc.py

"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

__author__ = 'annacuomo'


SAIGE_QTL_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/saige-qtl:8d840248fc025d5f58577ef03a6c23634ee941e5'


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
    saige_command_step3 = '/usr/local/bin/Rscript step3_gene_pvalue_qtl.R'
    saige_command_step3 += f' --assocFile={saige_sv_output_file}'
    saige_command_step3 += f' --geneName={gene_name}'
    saige_command_step3 += f' --genePval_outputFile={saige_gene_pval_output_file}'
    return saige_command_step3


config = get_config()


@click.command()
@click.option(
    '--pheno-cov-filename',
    default='/usr/local/bin/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt',
)
@click.option('--vcf-file-path', default='/usr/local/bin/genotype_10markers.vcf.gz')
@click.option('--vcf-field', default='GT')
@click.option('--covs-list', default='X1,X2,pf1,pf2')
@click.option('--sample-covs-list', default='X1,X2')
@click.option('--sample-id', default='IND_ID')
@click.option(
    '--null-output-path',
    default=output_path('nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1'),
)
@click.option(
    '--sv-output-path',
    default=output_path('nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis'),
)
@click.option(
    '--gene-pvals-output-path',
    default=output_path(
        'nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_genePval'
    ),
)
@click.option('--plink-path', default='/usr/local/bin/n.indep_100_n.cell_1_01.step1')
@click.option('--gene-name', default='gene_1')
@click.option('--chrom', default='2')
@click.option('--cis-window-file', default='/usr/local/bin/gene_1_cis_region.txt')
@click.option('--fit-null-mem', default='10Gi')
@click.option('--run-sv-assoc-mem', default='10Gi')
@click.option('--get-gene-pvals-mem', default='10Gi')
def association_pipeline(
    pheno_cov_filename: str,
    vcf_file_path: str,
    vcf_field: str,
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
    fit_null_mem: str,
    run_sv_assoc_mem: str,
    get_gene_pvals_mem: str,
):
    """
    Run association for one gene
    """

    batch = get_batch('SAIGE-QTL pipeline')

    # step 1 (fit null)
    fit_null_job = batch.new_job(name='fit null')
    fit_null_job.image(SAIGE_QTL_IMAGE)
    fit_null_job.storage(fit_null_mem)
    fit_null_job.declare_resource_group(
        output={
            'rda': '{root}.rda',
            'varianceRatio': '{root}.varianceRatio.txt',
        }
    )
    cmd = build_fit_null_command(
        pheno_file=pheno_cov_filename,
        cov_col_list=covs_list,
        sample_cov_col_list=sample_covs_list,
        sample_id_pheno=sample_id,
        output_prefix=fit_null_job.output,
        plink_path=plink_path,
        pheno_col=gene_name,
    )
    fit_null_job.command(cmd)

    # copy the output file to persistent storage?
    if null_output_path.startswith('gs://'):
        batch.write_output(fit_null_job.output, null_output_path)

    # step 2 (cis eQTL single variant test)
    run_sv_assoc_job = batch.new_job(name='single variant test')
    run_sv_assoc_job.image(SAIGE_QTL_IMAGE)
    run_sv_assoc_job.storage(run_sv_assoc_mem)
    run_sv_assoc_job.depends_on(fit_null_job)
    cmd = build_run_single_variant_test_command(
        vcf_file=vcf_file_path,
        vcf_file_index=f'{vcf_file_path}.csi',
        vcf_field=vcf_field,
        saige_output_file=run_sv_assoc_job.output,
        chrom=chrom,
        cis_window_file=cis_window_file,
        gmmat_model_path=fit_null_job.output['rda'],
        variance_ratio_path=fit_null_job.output['varianceRatio'],
    )
    run_sv_assoc_job.command(cmd)
    if sv_output_path.startswith('gs://'):
        batch.write_output(run_sv_assoc_job.output, sv_output_path)

    # step3 (gene-level pvalues)
    get_gene_pvals_job = batch.new_job(name='gene level pvalues')
    get_gene_pvals_job.image(SAIGE_QTL_IMAGE)
    get_gene_pvals_job.storage(get_gene_pvals_mem)
    get_gene_pvals_job.depends_on(run_sv_assoc_job)
    cmd = build_obtain_gene_level_pvals_command(
        gene_name=gene_name,
        saige_sv_output_file=run_sv_assoc_job.output,
        saige_gene_pval_output_file=get_gene_pvals_job.output,
    )
    get_gene_pvals_job.command(cmd)
    if gene_pvals_output_path.startswith('gs://'):
        batch.write_output(get_gene_pvals_job.output, gene_pvals_output_path)

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    association_pipeline()  # pylint: disable=no-value-for-parameter
