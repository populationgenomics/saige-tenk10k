#!/usr/bin/env python3
# pylint: disable=too-many-arguments

"""
Hail Batch workflow to perform association tests using SAIGE-QTL

- given all input files already generated
    - pheno_cov filename, genotypes with annotations
- builds saige commands (just add to str)
- run saige commands (execute Rscript from command line)
- aggregate & summarise results (not yet)

"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch

__author__ = 'annacuomo'


# this needs to be updated
SAIGE_QTL_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/saige-qtl'


# Fit null model (step 1)
def build_fit_null_command(
    pheno_file: str,
    cov_col_list: str,  # PC1,age,sex
    sample_cov_col_list: str,  # age,sex
    sample_id_pheno: str,  # individual_id
    output_prefix: str,  # ./${gene_name}_${celltype}_count
    plink_path: str,  # ./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly
    pheno_col: str = 'y',  # ${gene_name}, can be random y because we are doing one gene
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
    saige_command_step1 = 'Rscript step1_fitNULLGLMM_qtl.R'
    saige_command_step1 += f' --useSparseGRMtoFitNULL={use_sparse_grm_null}'
    saige_command_step1 += f' --useGRMtoFitNULL={use_grm_null}'
    saige_command_step1 += f' --phenoFile={pheno_file}'
    saige_command_step1 += f' --phenoCol={pheno_col}'
    saige_command_step1 += f' --covarColList={cov_col_list}'
    saige_command_step1 += f' --sampleCovarColList={sample_cov_col_list}'
    saige_command_step1 += f' --sampleIDColinphenoFile={sample_id_pheno}'
    saige_command_step1 += f' --traitType={trait_type}'
    saige_command_step1 += f' --outputPrefix={output_prefix}'
    saige_command_step1 += f' --skipVarianceRatioEstimation={skip_vre}'
    saige_command_step1 += f' --isRemoveZerosinPheno={pheno_remove_zeros}'
    saige_command_step1 += f' --isCovariateOffset={is_cov_offset}'
    saige_command_step1 += f' --isCovariateTransform={is_cov_transform}'
    saige_command_step1 += f' --skipModelFitting={skip_model_fitting}'
    saige_command_step1 += f' --tol={tol}'
    saige_command_step1 += f' --plinkFile={plink_path}'
    saige_command_step1 += f' --IsOverwriteVarianceRatioFile={is_overwrite_vre_file}'
    return saige_command_step1


# Run single variant association (step 2)
def build_run_single_variant_test_command(
    vcf_file: str,  # chr${i}.dose.filtered.R2_0.8.vcf.gz
    vcf_file_index: str,  # chr${i}.dose.filtered.R2_0.8.vcf.gz.csi
    saige_output_file: str,  # output/${prefix}_chr${i}${datatype}
    chrom: str,  # ${i}
    cis_window_file: str,
    gmmat_model_path: str,  # generated by step1 (${prefix}.rda)
    variance_ratio_path: str,  # generated by step1 (${prefix}.varianceRatio.txt)
    # allele_order: str = 'alt-first',  # change this and skip 2-g??
    min_maf: float = 0,
    min_mac: int = 5,
    loco_bool: str = 'FALSE',
    # is_fast_test: str = 'TRUE',
    n_markers: int = 10000,
    # pval_cutoff: float = 0.05,
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
    saige_command_step2_sv = 'Rscript step2_tests_qtl.R'
    saige_command_step2_sv += f' --vcfFile={vcf_file}'
    saige_command_step2_sv += f' --vcfFileIndex={vcf_file_index}'
    saige_command_step2_sv += f' --SAIGEOutputFile={saige_output_file}'
    saige_command_step2_sv += f' --chrom={chrom}'
    saige_command_step2_sv += f' --minMAF={min_maf}'
    saige_command_step2_sv += f' --minMAC={min_mac}'
    saige_command_step2_sv += f' --LOCO={loco_bool}'
    saige_command_step2_sv += f' --GMMATmodelFile={gmmat_model_path}'
    saige_command_step2_sv += f' --SPAcutoff={spa_cutoff}'
    saige_command_step2_sv += f' --varianceRatioFile={variance_ratio_path}'
    # saige_command_step2_sv += f' --is_fastTest={is_fast_test}'
    saige_command_step2_sv += f' --rangestoIncludeFile={cis_window_file}'
    saige_command_step2_sv += f' --markers_per_chunk={n_markers}'
    # saige_command_step2_sv += f' --pval_cutoff_for_fastTest={pval_cutoff}'
    return saige_command_step2_sv


def build_obtain_gene_level_pvals_command(
    gene_name: str,
    saige_sv_output_file: str,
    saige_gene_pval_output_file: str,
):
    """Build SAIGE command to obtain gene-level pvals
    Only for single-variant tests (Step3)
    combines single-variant p-values to obtain one gene
    level p-valye

    Input:
    - ouput of previous step, association file (txt)
    - gene we need to aggregate results for (across SNPs)
    - path for output file
    """
    saige_command_step3 = 'Rscript step3_gene_pvalue_qtl.R'
    saige_command_step3 += f' --assocFile={saige_sv_output_file}'
    saige_command_step3 += f' --geneName={gene_name}'
    saige_command_step3 += f' --genePval_outputFile={saige_gene_pval_output_file}'
    return saige_command_step3


# # Run gene set association (for RVs)
# def build_run_set_test_command(
#     vcf_file: str,  # chr${i}.dose.filtered.R2_0.8.vcf.gz
#     vcf_file_index: str,  # chr${i}.dose.filtered.R2_0.8.vcf.gz.tbi
#     saige_output_file: str,  # should end in txt?
#     chrom: str,  # check
#     gmmat_model_path: str,  # generated by step1 (.rda)
#     variance_ratio_path: str,  # generated by step1 (.txt)
#     group_annotation: str,  # e.g., 'lof:missense'
#     group_file: str,  # .txt
#     # allele_order: str = 'alt-first',  # change this and skip 2-g??
#     min_maf: float = 0,
#     min_mac: int = 5,
#     loco_bool: str = 'FALSE',
#     n_markers: int = 10000,
#     pval_cutoff: float = 0.05,
#     spa_cutoff: int = 2,
#     max_maf_group: str = '0.5',
#     min_maf_group: str = '0',
# ):
#     """Build SAIGE command for running set test
#     This will run a set test using Burden, SKAT and SKAT-O

#     Input:
#     plink_prefix: path to plink files (bim, bed, fam)
#     saige output path: path to output saige file
#     chrom: chromosome to run this on
#     GMMAT model file: null model fit from previous step (.rda)
#     Variance Ratio file: as estimated from previous step (.txt)
#     group annotation: select only specific annotations from group file (e.g., lof)
#     group file: for each gene/set,
#                 one row specifying variants,
#                 one row specifying each variant's anno,
#                 one optional row with all weights
#     allele order: specifying whether alt-first or ref-first in genotype files
#     min MAF: minimum variant minor allele frequency to include
#     min MAC: minimum variant minor allele count to include
#     LOCO: leave one chromosome out
#     is no adjusted cov: covariate adjustment?
#     specify whether we're using a sparse GRM (vs full? vs no GRM at all?)

#     Output:
#     Rscript command (str) ready to run
#     """
#     saige_command_step2 = 'Rscript step2_tests_qtl.R'
#     saige_command_step2 += f' --vcfFile={vcf_file}'
#     saige_command_step2 += f' --vcfFileIndex={vcf_file_index}'
#     # saige_command_step2 += f' --AlleleOrder={allele_order}'
#     saige_command_step2 += f' --SAIGEOutputFile={saige_output_file}'
#     saige_command_step2 += f' --chrom={chrom}'
#     saige_command_step2 += f' --minMAF={min_maf}'
#     saige_command_step2 += f' --minMAC={min_mac}'
#     saige_command_step2 += f' --LOCO={loco_bool}'
#     saige_command_step2 += f' --GMMATmodelFile={gmmat_model_path}'
#     saige_command_step2 += f' --SPAcutoff={spa_cutoff}'
#     saige_command_step2 += f' --varianceRatioFile={variance_ratio_path}'
#     saige_command_step2 += f' --groupFile={group_file}'
#     saige_command_step2 += f' --annotation_in_groupTest={group_annotation}'
#     saige_command_step2 += f' --maxMAF_in_groupTest={max_maf_group}'
#     saige_command_step2 += f' --minMAF_in_groupTest={min_maf_group}'
#     # saige_command_step2 += f' --is_noadjCov={is_no_adj_cov}'
#     # saige_command_step2 += f' --is_sparseGRM={is_sparse_grm}' # related individuals
#     saige_command_step2 += f' --markers_per_chunk={n_markers}'
#     saige_command_step2 += f' --pval_cutoff_for_fastTest={pval_cutoff}'
#     # saige_command_step2 += f' --is_EmpSPA={is_emp_spa}'

#     return saige_command_step2


config = get_config()


@click.command()
@click.option(
    '--pheno-cov-filename',
    default='input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt',
)
@click.option('--vcf-file-path', default='input/genotype_10markers.vcf.gz')
@click.option('--covs-list', default='X1,X2,pf1,pf2')
@click.option('--sample-covs-list', default='X1,X2')
@click.option('--sample-id', default='IND_ID')
@click.option(
    '--null-output-path',
    default='output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1',
)
@click.option(
    '--sv-output-path',
    default='output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis',
)
@click.option(
    '--gene-pvals-output-path',
    default='output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_genePval',
)
@click.option('--plink-path', default='input/n.indep_100_n.cell_1_01.step1')
@click.option('--gene-name', default='gene_1')
@click.option('--chrom', default='2')
@click.option('--cis-window-file', default='input/gene_1_cis_region.txt')
def association_pipeline(
    pheno_cov_filename_tsv: str,
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
    # chromosomes: str,
    # celltypes: str,
):
    """
    Run association for one gene
    """

    batch = get_batch('SAIGE-QTL pipeline')

    # step 1 (fit null)
    fit_null_job = batch.new_job(name='fit null')
    fit_null_job.image(SAIGE_QTL_IMAGE)
    cmd = build_fit_null_command(
        pheno_file=pheno_cov_filename_tsv,
        cov_col_list=covs_list,
        sample_cov_col_list=sample_covs_list,
        sample_id_pheno=sample_id,
        output_prefix=null_output_path,
        plink_path=plink_path,
        pheno_col=gene_name,
    )
    fit_null_job.command(cmd)

    # step 2 (cis eQTL single variant test)
    run_sv_assoc_job = batch.new_job(name='single variant test')
    run_sv_assoc_job.image(SAIGE_QTL_IMAGE)
    run_sv_assoc_job.depends_on(fit_null_job)
    cmd = build_run_single_variant_test_command(
        vcf_file=vcf_file_path,
        vcf_file_index=f'{vcf_file_path}.csi',
        saige_output_file=sv_output_path,
        chrom=chrom,
        cis_window_file=cis_window_file,
        gmmat_model_path=f'{null_output_path}.rda',
        variance_ratio_path=f'{null_output_path}.varianceRatio.txt',
    )
    run_sv_assoc_job.command(cmd)

    # step3 (gene-level pvalues)
    get_gene_pvals_job = batch.new_job(name='gene level pvalues')
    get_gene_pvals_job.image(SAIGE_QTL_IMAGE)
    get_gene_pvals_job.depends_on(run_sv_assoc_job)
    cmd = build_obtain_gene_level_pvals_command(
        gene_name=gene_name,
        saige_sv_output_file=sv_output_path,
        saige_gene_pval_output_file=gene_pvals_output_path,
    )
    get_gene_pvals_job.command(cmd)

    # celltypes_list = celltypes.split(' ')
    # chromosomes_list = chromosomes.split(' ')

    # for celltype in celltypes_list:
    #     for chromosome in chromosomes_list:
    #         # get genes
    #         for gene in genes:

    #             fit_null_job = batch.new_job(name='fit null')
    #             fit_null_job.image(SAIGE_QTL_IMAGE)
    #             cmd = build_fit_null_command(
    #                 pheno_file=pheno_cov_filename_tsv,
    #                 cov_col_list=covs_list,
    #                 sample_id_pheno=sample_id,
    #                 output_prefix=output_path,
    #                 plink_path=plink_path,
    #                 pheno_col=pheno_col,
    #             )
    #             fit_null_job.command(cmd)

    #             run_sv_assoc_job = batch.new_job(name='single variant test')
    #             run_sv_assoc_job.image(SAIGE_QTL_IMAGE)
    #             run_sv_assoc_job.depends_on(fit_null_job)
    #             cmd = build_run_single_variant_test_command(
    #                 vcf_file=vcf_file_path,
    #                 vcf_file_index=f'{vcf_file_path}.tbi',
    #                 saige_output_file=output_path,
    #             )
    #             run_sv_assoc_job.command(cmd)

    #             run_set_assoc_job = batch.new_job(name='settest')
    #             run_set_assoc_job.image(SAIGE_QTL_IMAGE)
    #             run_set_assoc_job.depends_on(fit_null_job)
    #             cmd = build_run_set_test_command(
    #                 vcf_file=vcf_file_path,
    #                 vcf_file_index=f'{vcf_file_path}.tbi',
    #                 saige_output_file=output_path,
    #             )
    #             run_set_assoc_job.command(cmd)

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    association_pipeline()  # pylint: disable=no-value-for-parameter
