#!/usr/bin/env python3
# pylint: disable=import-error

__author__ = 'annacuomo'

"""
Hail Batch workflow to perform association tests using SAIGE-QTL

- given all input files already generated (pheno_cov filename, genotypes with annotations)
- builds and runs saige commands
- summarise results

"""

# Fit null model
def build_fit_null_command(
    sparse_grm_file: str,  # data/input/nfam_5_nindep_0.mtx
    sparse_grm_sampleid_file: str,  # data/input/nfam_5_nindep_0.mtx.sampleID
    pheno_file: str,
    cov_col_list: str,  # PC1
    sample_id_pheno: str,  # IND_ID
    plink_path: str,  # ./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly
    output_prefix: str,
    pheno_col: str = 'y',
    trait_type: str = 'count',
    skip_vre: str = 'FALSE',  # this is a boolean but that's encoded differently between R and python
    is_cov_offset: str = 'FALSE',
    is_cov_transform: str = 'FALSE',
    skip_model_fitting: str = 'FALSE',
    is_overwrite_vre_file: str = 'TRUE',
):
    """Build SAIGE command for fitting null model
    This will fit a Poisson / NB mixed model under the null hypothesis

    Input:
    - Sparse GRM file (generated previously using build_sparse_grm_command)
    - Sparse GRM sample ID file (to match sample ID's? check)
    - Phenotype / covariate file - rows: samples, cols: pheno (y), cov1, cov2 etc
    - List of columns from previous file that should be used as covariates
    - Column name of sample column (e.g., IND_ID, or CPG_ID etc)
    - Plink path: path to plink file
    - output prefix: where to save the fitted model (.rda)
    - pheno col: name of column specifying pheno (default: "y")
    - trait type: count = Poisson, count_nb = Negative Binomial
    - option to skip Variance Ratio estimation (discouraged)
    - option to add an offset to the fixed covariates (or if there are no covariates add an intercept?)
    - option to transform (scale?) covariates?
    - option to skip model fitting (discouraged)
    - genotype file to estimate variance ratio (plink format) - which variants?
    - overwrite variance ratio file (estimated here)

    Output:
    Rscript command (str) ready to run (bash)
    """
    saige_command_step1 = 'Rscript step1_fitNULLGLMM_qtl.R'
    # figure out whether these are always needed or no
    saige_command_step1 += f' --sparseGRMFile={sparse_grm_file}'
    saige_command_step1 += f' --sparseGRMSampleIDFile={sparse_grm_sampleid_file}'
    saige_command_step1 += f' --phenoFile={pheno_file}'
    saige_command_step1 += f' --phenoCol={pheno_col}'
    saige_command_step1 += f' --covarColList={cov_col_list}'
    saige_command_step1 += f' --sampleIDColinphenoFile={sample_id_pheno}'
    saige_command_step1 += f' --traitType={trait_type}'
    saige_command_step1 += f' --outputPrefix={output_prefix}'
    saige_command_step1 += f' --skipVarianceRatioEstimation={skip_vre}'
    saige_command_step1 += f' --isCovariateOffset={is_cov_offset}'
    saige_command_step1 += f' --isCovariateTransform={is_cov_transform}'
    saige_command_step1 += f' --skipModelFitting={skip_model_fitting}'
    saige_command_step1 += f' --plinkFile={plink_path}'
    saige_command_step1 += f' --IsOverwriteVarianceRatioFile={is_overwrite_vre_file}'
    return saige_command_step1


# Run gene set association
def build_run_set_test_command(
    plink_prefix: str,  # ./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly
    saige_output_file: str,  # should end in txt?
    chrom: str,  # check
    gmmat_model_path: str,  # generated by step1 (.rda)
    variance_ratio_path: str,  # generated by step1 (.txt)
    group_annotation: str,  # e.g., 'lof:missense'
    group_file: str,  # .txt
    allele_order: str = 'alt-first',  # change this and skip 2-g??
    min_maf: float = 0,
    min_mac: int = 5,
    loco_bool: str = 'FALSE',
    is_no_adj_cov: str = 'TRUE',
    is_sparse_grm: str = 'FALSE',
    n_markers: int = 10000,
    pval_cutoff: float = 0.05,
    spa_cutoff: int = 10000,
    is_emp_spa: str = 'FALSE',
    max_maf_group: float = 0.5,
):
    """Build SAIGE command for running set test
    This will run a set test using Burden, SKAT and SKAT-O

    Input:
    plink_prefix: path to plink files (bim, bed, fam)
    saige ouput path: path to output saige file
    chrom: chromosome to run this on
    GMMAT model file: null model fit from previous step (.rda)
    Variance Ratio file: as estimated from previous step (.txt)
    group annotation: select only specific annotations from group file (e.g., lof)
    group file: for each gene/set, one row specifying variants, one row specifying each variant's anno, one optional row with all weights
    allele order: specifying whether alt-first or ref-first in genotype files
    min MAF: minimum variant minor allele frequency to include
    min MAC: minimum variant minor allele count to include
    LOCO: leave one chromosome out (for what specifically?)
    is no adjusted cov: covariate adjustment?
    specify whether we're using a sparse GRM (vs full? vs no GRM at all?)

    Output:
    Rscript command (str) ready to run
    """
    saige_command_step2 = 'Rscript step2_tests_qtl.R'
    saige_command_step2 += f' --bedFile={plink_prefix}.bed'
    saige_command_step2 += f' --bimFile={plink_prefix}.bim'
    saige_command_step2 += f' --famFile={plink_prefix}.fam'
    saige_command_step2 += f' --AlleleOrder={allele_order}'
    saige_command_step2 += f' --SAIGEOutputFile={saige_output_file}'
    saige_command_step2 += f' --chrom={chrom}'
    saige_command_step2 += f' --minMAF={min_maf}'
    saige_command_step2 += f' --minMAC={min_mac}'
    saige_command_step2 += f' --LOCO={loco_bool}'
    saige_command_step2 += f' --GMMATmodelFile={gmmat_model_path}'
    saige_command_step2 += f' --varianceRatioFile={variance_ratio_path}'
    saige_command_step2 += f' --is_noadjCov={is_no_adj_cov}'
    saige_command_step2 += f' --is_sparseGRM={is_sparse_grm}'
    saige_command_step2 += f' --markers_per_chunk={n_markers}'
    saige_command_step2 += f' --pval_cutoff_for_fastTest={pval_cutoff}'
    saige_command_step2 += f' --SPAcutoff={spa_cutoff}'
    saige_command_step2 += f' --is_EmpSPA={is_emp_spa}'
    saige_command_step2 += f' --annotation_in_groupTest={group_annotation}'
    saige_command_step2 += f' --maxMAF_in_groupTest={max_maf_group}'
    saige_command_step2 += f' --groupFile={group_file}'
    return saige_command_step2