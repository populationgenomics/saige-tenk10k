#!/usr/bin/env python3
# pylint: disable=import-error,import-outside-toplevel,missing-function-docstring,no-value-for-parameter,too-many-arguments,too-many-locals,wrong-import-order,wrong-import-position,consider-using-enumerate,chained-comparison

__author__ = "annacuomo"

"""
Hail Batch workflow for the rare-variant association analysis, including:

- get relevant variants around a gene and export genotypes as plink files
- generate input files for association tests
- run association tests
"""
# now copied from https://github.com/populationgenomics/cellregmap-pipeline/edit/main/batch.py

# import python modules
import os
import re
import sys

import click
import logging
from typing import Dict

from google.cloud import storage

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    copy_common_env,
    dataset_path,
    get_config,
    init_batch,
    output_path,
    remote_tmpdir,
)

import numpy as np
import pandas as pd
import xarray as xr

from limix.qc import quantile_gaussianize
from scipy.stats import shapiro

import hail as hl
import hailtop.batch as hb


# use logging to print statements, display at info level
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stderr,
)

DEFAULT_JOINT_CALL_MT = dataset_path("mt/v7.mt")
DEFAULT_ANNOTATION_HT = dataset_path(
    "tob_wgs_vep/104/vep104.3_GRCh38.ht"
)  # atm VEP only


CELLREGMAP_IMAGE = get_config()["workflow"][
    "driver_image"
]  # australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:dev
SAIGE_QTL_IMAGE = "australia-southeast1-docker.pkg.dev/cpg-common/images/saige-qtl"  # do I need the get_config part?

MULTIPY_IMAGE = "australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16"  # not sure I will need this

# region SUBSET_VARIANTS

# same as https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py
def filter_variants(
    mt_path: str,  # "mt/v7.mt"
    samples: list[str],
    output_mt_path: str,  # "tob_wgs_rv/densified_rv_only.mt"
    grm_plink_file: str,
):
    """Subset hail matrix table

    Input:
    joint call hail matrix table
    set of samples for which we have scRNA-seq data

    Output:
    subset hail matrix table, containing only variants that:
    1) are not ref-only, 2) biallelic, 3) meet QC filters, 4) are rare (MAF<5%)

    also, plink file containing variants that satisfy 1),2),3)
    but that are common (MAF>1%) to build sparse GRM
    """
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # subset to relevant samples (samples we have scRNA-seq data for)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # densify
    mt = hl.experimental.densify(mt)

    # filter out low quality variants and consider biallelic SNPss only (no multi-allelic, no ref-only, no indels)
    mt = mt.filter_rows(  # check these filters!
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs
    )

    mt = hl.variant_qc(mt)
    # filter common (enough) variants to build sparse GRM
    grm_mat = mt.filter_rows(
        (mt.variant_qc.AF[1] > 0.01) & (mt.variant_qc.AF[1] < 1)
        | (mt.variant_qc.AF[1] < 0.99) & (mt.variant_qc.AF[1] > 0)
    )

    # export to plink
    from hail.methods import export_plink

    export_plink(grm_mt, grm_plink_file, ind_id=grm_mt.s)

    # filter rare variants only (MAF < 5%)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.AF[1] > 0)
        | (mt.variant_qc.AF[1] > 0.95) & (mt.variant_qc.AF[1] < 1)
    )
    mt.write(output_mt_path, overwrite=True)
    logging.info(f"Number of rare (freq<5%) and QCed biallelic SNPs: {mt.count()[0]}")


# endregion SUBSET_VARIANTS

# region CREATE_SPARSE_GRM

# NEW
def build_sparse_grm_command(
    plink_path: str,  # ./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly
    output_prefix: str,  # should end in sparseGRM
    n_threads: int = 4,
    n_random: int = 2000,
    relatedness_cutoff: float = 0.125,
):
    """Build SAIGE command for SPARSE GRM

    Input:
    plink_path: path to (GRM) plink file
    various flags

    Output:
    Rscript command (str) ready to run
    """
    saige_command_step0 = "Rscript createSparseGRM.R"
    saige_command_step0 += f" --plinkFile={plink_path}"
    saige_command_step0 += f" --nThreads={n_threads}"
    saige_command_step0 += f" --outputPrefix={output_prefix}"
    saige_command_step0 += f" --numRandomMarkerforSparseKin={n_random}"
    saige_command_step0 += f" --relatednessCutoff={relatedness_cutoff}"
    return saige_command_step0


# endregion CREATE_SPARSE_GRM

# region GET_GENE_SPECIFIC_VARIANTS

# same as https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py
def get_promoter_variants(
    mt_path: str,  # ouput path from function above
    ht_path: str,
    gene_details: dict[str, str],  # output of make_gene_loc_dict
    window_size: int,
    plink_file: str,  # "tob_wgs_rv/pseudobulk_rv_association/plink_files/GENE"
):
    """Subset hail matrix table

    Input:
    mt_path: path to already subsetted hail matrix table
    ht_path: path to VEP HT
    gene_details: dict of info for current gene
    window_size: int, size of flanking region around genes
    plink_file: str, file prefix for writing plink data

    Output:
    For retained variants, that are: 1) in promoter regions and
    2) within 50kb up or down-stream of the gene body (or in the gene body itself)
    (on top of all filters done above)

    returns nothing
    """

    # read hail matrix table object (pre-filtered)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    gene_name = gene_details["gene_name"]

    # get relevant chromosome
    chrom = gene_details["chr"]

    # subset to window
    # get gene body position (start and end) and build interval
    left_boundary = max(1, int(gene_details["start"]) - window_size)
    right_boundary = min(
        int(gene_details["end"]) + window_size,
        hl.get_reference("GRCh38").lengths[chrom],
    )
    # get gene-specific genomic interval
    gene_interval = f"{chrom}:{left_boundary}-{right_boundary}"
    logging.info(f"Interval considered: {gene_interval}")  # "chr22:23219960-23348287"

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome="GRCh38")]
    )
    mt_path = output_path(f"{gene_name}_in_window.mt", "tmp")
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f"Number of variants within interval: {mt.count()[0]}")

    # annotate using VEP
    vep_ht = hl.read_table(ht_path)
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)

    # filter variants found to be in promoter regions
    mt = mt.filter_rows(
        mt.vep.regulatory_feature_consequences["biotype"].contains("promoter")
    )
    promoter_path = output_path(f"{gene_name}promoter_variants.mt", "tmp")
    mt = mt.checkpoint(promoter_path, overwrite=True)  # checkpoint
    logging.info(
        f"Number of rare (freq<5%) QC-passing, biallelic SNPs in promoter regions: {mt.count()[0]}"
    )

    # export this as a Hail table for downstream analysis
    # consider exporting the whole MT?
    ht_path = output_path(
        f"summary_hts/{gene_name}_rare_promoter_summary.ht", "analysis"
    )
    ht = mt.rows()
    ht.write(ht_path, overwrite=True)

    # export MT object to PLINK (promoter variants)
    # pylint: disable=import-outside-toplevel
    from hail.methods import export_plink

    export_plink(mt, plink_file, ind_id=mt.s)


# endregion GET_GENE_SPECIFIC_VARIANTS

# region PREPARE_INPUT_FILES

# this will need to change:
# pheno and cov may be included within the same file
# geno can be plink directly
# kinship is not necessary, can be computed using step0 of saige (create sparse grm)
# a group file will be needed (for gene-set tests)
def prepare_input_files(
    gene_name: str,
    cell_type: str,
    genotype_file_bed: str,
    genotype_file_bim: str,
    genotype_file_fam: str,
    phenotype_file: str,
    kinship_file: str,
    sample_mapping_file: str,
):
    """Prepare association test input files

    Input:
    genotype: plink files
    phenotype: gene expression

    Output:
    genotype matrix
    phenotype vector
    """

    # check that variants exist for this gene, before importing
    bim_file = to_path(genotype_file_bim)
    if (not bim_file.exists()) or bim_file.stat().st_size == 0:
        logging.info(f"{bim_file} does not exist")
        return None

    from pandas_plink import read_plink1_bin

    expression_filename = to_path(
        output_path(f"expression_files/{gene_name}_{cell_type}.csv")
    )
    genotype_filename = to_path(
        output_path(f"genotype_files/{gene_name}_rare_regulatory.csv")
    )
    kinship_filename = to_path(output_path(f"{gene_name}_kinship_common_samples.csv"))

    # read in phenotype file (tsv)
    phenotype = pd.read_csv(phenotype_file, sep="\t", index_col=0)

    phenotype = xr.DataArray(
        phenotype.values,
        dims=["sample", "gene"],
        coords={"sample": phenotype.index.values, "gene": phenotype.columns.values},
    )

    # read in genotype file (plink format)
    to_path(genotype_file_bed).copy("temp.bed")  # bed
    bim_file.copy("temp.bim")  # bim
    to_path(genotype_file_fam).copy("temp.fam")  # fam

    geno = read_plink1_bin("temp.bed")

    if kinship_file is not None:
        # read in GRM (genotype relationship matrix; kinship matrix)
        kinship = pd.read_csv(kinship_file, index_col=0)
        kinship.index = kinship.index.astype("str")
        assert all(
            kinship.columns == kinship.index
        )  # symmetric matrix, donors x donors
        kinship = xr.DataArray(
            kinship.values,
            dims=["sample_0", "sample_1"],
            coords={"sample_0": kinship.columns, "sample_1": kinship.index},
        )
        kinship = kinship.sortby("sample_0").sortby("sample_1")

    # this file will map different IDs (and OneK1K ID to CPG ID)
    sample_mapping = pd.read_csv(dataset_path(sample_mapping_file), sep="\t")

    # ensure samples are the same and in the same order across input files
    # samples with expression data
    donors_exprs = set(phenotype.sample.values).intersection(
        set(sample_mapping["OneK1K_ID"].unique())
    )

    logging.info(f"Number of unique donors with expression data: {len(donors_exprs)}")

    # samples with genotype data
    donors_geno = set(geno.sample.values).intersection(
        set(sample_mapping["InternalID"].unique())
    )
    logging.info(f"Number of unique donors with genotype data: {len(donors_geno)}")

    # samples with both (can this be done in one step?)
    sample_mapping1 = sample_mapping.loc[sample_mapping["OneK1K_ID"].isin(donors_exprs)]
    sample_mapping_both = sample_mapping1.loc[
        sample_mapping1["InternalID"].isin(donors_geno)
    ]
    donors_e = sample_mapping_both["OneK1K_ID"].unique()
    donors_g = sample_mapping_both["InternalID"].unique()
    assert len(donors_e) == len(donors_g)

    if kinship_file is not None:
        # samples in kinship
        donors_e_short = [re.sub(".*_", "", donor) for donor in donors_e]
        donors_k = sorted(
            set(list(kinship.sample_0.values)).intersection(donors_e_short)
        )

    logging.info(f"Number of unique common donors: {len(donors_g)}")

    # subset files

    # phenotype
    phenotype = phenotype.sel(sample=donors_e)
    # select gene
    y = phenotype.sel(gene=gene_name)
    y = quantile_gaussianize(y)
    del phenotype  # delete to free up memory
    # make data frame to save as csv
    y_df = pd.DataFrame(
        data=y.values.reshape(y.shape[0], 1), index=y.sample.values, columns=[gene_name]
    )

    # genotype
    geno = geno.sel(sample=donors_g)
    # make data frame to save as csv
    data = geno.values
    geno_df = pd.DataFrame(data, columns=geno.snp.values, index=geno.sample.values)
    geno_df = geno_df.dropna(axis=1)
    # delete large files to free up memory
    del geno

    if kinship_file is not None:
        # kinship
        kinship = kinship.sel(sample_0=donors_k, sample_1=donors_k)
        assert all(kinship.sample_0 == donors_k)
        assert all(kinship.sample_1 == donors_k)
        # make data frame to save as csv
        kinship_df = pd.DataFrame(
            kinship.values, columns=kinship.sample_0, index=kinship.sample_1
        )
        del kinship  # delete kinship to free up memory

    # save files
    with expression_filename.open("w") as ef:
        y_df.to_csv(ef, index=False)

    with genotype_filename.open("w") as gf:
        geno_df.to_csv(gf, index=False)

    if kinship_file is not None:
        with kinship_filename.open("w") as kf:
            kinship_df.to_csv(kf, index=False)
    else:
        kinship_df = None

    return y_df, geno_df, kinship_df


# endregion PREPARE_INPUT_FILES


# region GET_SAIGE_COMMANDS

# NEW
def build_fit_null_command(
    sparse_grm_file: str,  # data/input/nfam_5_nindep_0.mtx
    sparse_grm_sampleid_file: str,  # data/input/nfam_5_nindep_0.mtx.sampleID
    pheno_file: str,
    cov_col_list: str,
    sample_id_pheno: str,  # IND_ID
    plink_path: str,  # ./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly
    output_prefix: str,
    pheno_col: str = "y",
    trait_type: str = "count",
    skip_vre: str = "FALSE",  # this is a boolean but that's encoded differently between R and python
    is_cov_offset: str = "FALSE",
    is_cov_transform: str = "FALSE",
    skip_model_fitting: str = "FALSE",
    is_overwrite_vre_file: str = "TRUE",
):
    """Build SAIGE command for SPARSE GRM

    Input:
    - Sparse GRM file (generated previously using build_sparse_grm_command)
    - Sparse GRM sample ID file (to match sample ID's? check)
    - Phenotype / covariate file - rows: samples, cols: pheno (y), cov1, cov2 etc
    - List of columns from previous file that should be used as covariates
    - Column name of sample column (e.g., IND_ID, or CPG_ID etc)
    - Plink path: path to plink file
    - output prefix: where to save the fitted model (.rda)
    - pheno col: name of column specifying pheno (default: "y")

    Output:
    Rscript command (str) ready to run (bash)
    """
    saige_command_step1 = "Rscript step1_fitNULLGLMM_qtl.R"
    saige_command_step1 += f" --sparseGRMFile={sparse_grm_file}"
    saige_command_step1 += f" --sparseGRMSampleIDFile={sparse_grm_sampleid_file}"
    saige_command_step1 += f" --phenoFile={pheno_file}"
    saige_command_step1 += f" --phenoCol={pheno_col}"
    saige_command_step1 += f" --covarColList={cov_col_list}"
    saige_command_step1 += f" --sampleIDColinphenoFile={sample_id_pheno}"
    saige_command_step1 += f" --traitType={trait_type}"
    saige_command_step1 += f" --outputPrefix={output_prefix}"
    saige_command_step1 += f" --skipVarianceRatioEstimation={skip_vre}"
    saige_command_step1 += f" --isCovariateOffset={is_cov_offset}"
    saige_command_step1 += f" --isCovariateTransform={is_cov_transform}"
    saige_command_step1 += f" --skipModelFitting={skip_model_fitting}"
    saige_command_step1 += f" --plinkFile={plink_path}"
    saige_command_step1 += f" --IsOverwriteVarianceRatioFile={is_overwrite_vre_file}"
    return saige_command_step1

# NEW
def build_run_set_test_command(
    plink_prefix: str,  # ./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly
    saige_output_file: str,  # should end in txt?
    chrom: str,  # check
    gmmat_model_path: str,  # generated by step1 (.rda)
    variance_ratio_path: str,  # generated by step1 (.txt)
    group_annotation: str,  # e.g., "lof:missense"
    group_file: str,  # .txt
    allele_order: str = "alt-first",  # change this and skip 2-g??
    min_maf: float = 0,
    min_mac: int = 5,
    loco_bool: str = "FALSE",
    is_no_adj_cov: str = "TRUE",
    is_sparse_grm: str = "FALSE",
    n_markers: int = 10000,
    pval_cutoff: float = 0.05,
    spa_cutoff: int = 10000,
    is_emp_spa: str = "FALSE",
    max_maf_group: float = 0.5,
):
    """Build SAIGE command for SPARSE GRM

    Input:
    plink_path: path to (GRM) plink file
    various flags

    Output:
    Rscript command (str) ready to run
    """
    saige_command_step2 = "Rscript step2_tests_qtl.R"
    saige_command_step2 += f" --bedFile={plink_prefix}.bed"
    saige_command_step2 += f" --bimFile={plink_prefix}.bim"
    saige_command_step2 += f" --famFile={plink_prefix}.fam"
    saige_command_step2 += f" --AlleleOrder={allele_order}"
    saige_command_step2 += f" --SAIGEOutputFile={saige_output_file}"
    saige_command_step2 += f" --chrom={chrom}"
    saige_command_step2 += f" --minMAF={min_maf}"
    saige_command_step2 += f" --minMAC={min_mac}"
    saige_command_step2 += f" --LOCO={loco_bool}"
    saige_command_step2 += f" --GMMATmodelFile={gmmat_model_path}"
    saige_command_step2 += f" --varianceRatioFile={variance_ratio_path}"
    saige_command_step2 += f" --is_noadjCov={is_no_adj_cov}"
    saige_command_step2 += f" --is_sparseGRM={is_sparse_grm}"
    saige_command_step2 += f" --markers_per_chunk={n_markers}"
    saige_command_step2 += f" --pval_cutoff_for_fastTest={pval_cutoff}"
    saige_command_step2 += f" --SPAcutoff={spa_cutoff}"
    saige_command_step2 += f" --is_EmpSPA={is_emp_spa}"
    saige_command_step2 += f" --annotation_in_groupTest={group_annotation}"
    saige_command_step2 += f" --maxMAF_in_groupTest={max_maf_group}"
    saige_command_step2 += f" --groupFile={group_file}"
    return saige_command_step2


# endregion GET_SAIGE_COMMANDS

# region RUN_ASSOCIATION


def run_gene_association(
    gene_name: str,  # "VPREB3"
    prepared_inputs: hb.resource.PythonResult,
    pv_path: str,  # "Bnaive/VPREB3_pvals.csv"
):
    """Run gene-set association test

    change
    """

    # if the previous method depdency returned None, we will fail to
    # unpack dataframes from it
    if prepared_inputs is None:
        return prepared_inputs

    from numpy import eye, ones

    # read the 3 dataframes generated by the previous job
    p_df, g_df, _ = prepared_inputs

    # because the current matrix is counting the copies of the reference allele
    # while we are interested in the alternative allele, flip the genotypes
    genotypes = 2 - g_df

    # get phenotypes
    pheno = p_df.values

    # covariates (intercept at present)
    covs = ones((genotypes.shape[0], 1))  # intercept of ones as covariates

    # contexts (no context-specific analysis now, just identity)
    contexts = eye(genotypes.shape[0])

    # create p-values data frame
    pvalues = get_crm_pvs(pheno, covs, genotypes, contexts)
    pv_df = pd.DataFrame(
        data=pvalues.reshape(pvalues.shape[0], 1).T,
        columns=cols,
        index=[gene_name],
    )

    pv_filename = to_path(pv_path)
    with pv_filename.open("w") as pf:
        pv_df.to_csv(pf)

    return str(pv_filename)


# endregion RUN_ASSOCIATION


# region AGGREGATE_RESULTS


def summarise_association_results(
    # pv_dfs: list[str],
    celltype: str,
    pv_all_filename_str: str,
):
    """Summarise results

    Input:
    p-values from all association tests

    Output:
    one csv table per cell type,
    combining results across all genes in a single file
    """
    from multipy.fdr import qvalue

    logging.info("before glob (pv files) - summarise job")
    storage_client = storage.Client()
    bucket = get_config()["storage"]["default"]["default"].removeprefix("gs://")
    prefix = f"{get_config()["workflow"]["output_prefix"]}/{celltype}/"
    existing_pv_files = set(
        f"gs://{bucket}/{filepath.name}"
        for filepath in storage_client.list_blobs(bucket, prefix=prefix, delimiter="/")
        if filepath.name.endswith("_pvals.csv")
    )
    logging.info(f"after glob - {len(existing_pv_files)} pv files to summarise")

    if len(existing_pv_files) == 0:
        raise Exception("No PV files, nothing to do")

    pv_all_df = pd.concat(
        [pd.read_csv(to_path(pv_df), index_col=0) for pv_df in existing_pv_files]
    )

    # run qvalues for all tests (multiple testing correction)
    _, qvals = qvalue(pv_all_df["P_CRM_VC"])
    pv_all_df["Q_CRM_VC"] = list(qvals)
    _, qvals = qvalue(pv_all_df["P_CRM_burden_max"])
    pv_all_df["Q_CRM_burden_max"] = list(qvals)
    _, qvals = qvalue(pv_all_df["P_CRM_burden_sum"])
    pv_all_df["Q_CRM_burden_sum"] = list(qvals)
    _, qvals = qvalue(pv_all_df["P_CRM_burden_comphet"])
    pv_all_df["Q_CRM_burden_comphet"] = list(qvals)
    _, qvals = qvalue(pv_all_df["P_CRM_omnibus_max"])
    pv_all_df["Q_CRM_omnibus_max"] = list(qvals)
    _, qvals = qvalue(pv_all_df["P_CRM_omnibus_sum"])
    pv_all_df["Q_CRM_omnibus_sum"] = list(qvals)
    _, qvals = qvalue(pv_all_df["P_CRM_omnibus_comphet"])
    pv_all_df["Q_CRM_omnibus_comphet"] = list(qvals)

    pv_all_filename = to_path(pv_all_filename_str)
    logging.info(f"Write summary results to {pv_all_filename}")
    with pv_all_filename.open("w") as pf:
        pv_all_df.to_csv(pf)


# endregion AGGREGATE_RESULTS

# region MISCELLANEOUS


def make_gene_loc_dict(file) -> dict[str, dict]:
    """
    Turn gene information into a dictionary
    to avoid opening this file for every gene
    """
    from csv import DictReader

    gene_dict = {}

    with open(to_path(file)) as handle:
        reader = DictReader(handle, delimiter="\t")

        for row in reader:
            gene_dict[row["gene_name"]] = row

    return gene_dict


# can probably be merged with below
def extract_genes(gene_list, expression_tsv_path) -> list[str]:
    """
    Takes a list of all genes and subsets to only those
    present in the expression file of interest
    """
    expression_df = pd.read_csv(to_path(expression_tsv_path), sep="\t")
    expression_df = filter_lowly_expressed_genes(expression_df)
    gene_ids = set(list(expression_df.columns.values)[1:])
    genes = set(gene_list).intersection(gene_ids)

    logging.info(f"Total genes to run: {len(list(sorted(genes)))}")

    return list(sorted(genes))


# copied from https://github.com/populationgenomics/tob-wgs/blob/main/scripts/eqtl_hail_batch/launch_eqtl_spearman.py
# generalised to specify min pct samples as input
def filter_lowly_expressed_genes(expression_df, min_pct=10):
    """Remove genes with low expression in all samples
    Input:
    expression_df: a data frame with samples as rows and genes as columns,
    containing normalised expression values (i.e., the average number of molecules
    for each gene detected in each person).
    Returns:
    A filtered version of the input data frame, after removing columns (genes)
    with 0 values in more than {min_pct}% of the rows (samples) - 10 by default.
    """
    # Remove genes with 0 expression in all samples
    expression_df = expression_df.loc[:, (expression_df != 0).any(axis=0)]
    genes_not_equal_zero = expression_df.iloc[:, 1:].values != 0
    n_expr_over_zero = pd.DataFrame(genes_not_equal_zero.sum(axis=0))
    percent_expr_over_zero = (n_expr_over_zero / len(expression_df.index)) * 100
    percent_expr_over_zero.index = expression_df.columns[1:]

    # Filter genes with less than 10 percent individuals with non-zero expression
    atleastNpercent = percent_expr_over_zero[(percent_expr_over_zero > min_pct)[0]]
    sample_ids = expression_df["sampleid"]
    expression_df = expression_df[atleastNpercent.index]
    expression_df.insert(loc=0, column="sampleid", value=sample_ids)

    return expression_df


def remove_sc_outliers(df, outliers=None):
    """
    Remove outlier samples, as identified by single-cell analysis
    """
    if outliers is None:
        outliers = ["966_967", "88_88"]
    else:
        outliers = outliers.extend(["966_967", "88_88"])
    df = df[-df["OneK1K_ID"].isin(outliers)]

    return df


# endregion MISCELLANEOUS

config = get_config()


@click.command()
@click.option("--celltypes")
@click.option("--expression-files-prefix", default="scrna-seq/grch38_association_files")
@click.option(
    "--sample-mapping-file-tsv",
    default="scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv",
)
@click.option("--mt-path", default=DEFAULT_JOINT_CALL_MT)
@click.option("--anno-ht-path", default=DEFAULT_ANNOTATION_HT)
@click.option(
    "--chromosomes",
    help="List of chromosome numbers to run rare variant association analysis on. "
    "Space separated, as one argument (Default: all)",
)
@click.option("--genes", default=None)
@click.option(
    "--max-gene-concurrency",
    type=int,
    default=50,
    help=(
        "To avoid resource starvation, set this concurrency to limit horizontal scale. "
        "Higher numbers have a better walltime, but risk jobs that are stuck (which are expensive)"
    ),
)
def saige_pipeline(
    celltypes: str,
    expression_files_prefix: str,
    sample_mapping_file_tsv: str,
    mt_path: str,
    anno_ht_path: str,
    chromosomes: str = "all",
    genes: str | None = None,
    window_size: int = 50000,
    max_gene_concurrency=100,
):

    sb = hb.ServiceBackend(
        billing_project=get_config()["hail"]["billing_project"],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch("CellRegMap pipeline", backend=sb)

    # extract samples for which we have single-cell (sc) data
    sample_mapping_file = pd.read_csv(dataset_path(sample_mapping_file_tsv), sep="\t")
    sample_mapping_file = remove_sc_outliers(sample_mapping_file)
    sc_samples = sample_mapping_file["InternalID"].unique()

    # filter to QC-passing, rare, biallelic variants
    output_mt_path = output_path("densified_rv_only.mt")
    if not to_path(output_mt_path).exists():

        filter_job = batch.new_python_job(name="MT filter job")
        copy_common_env(filter_job)
        filter_job.image(CELLREGMAP_IMAGE)
        filter_job.call(
            filter_variants,
            mt_path=mt_path,
            samples=list(sc_samples),
            output_mt_path=output_mt_path,
        )

    else:
        logging.info("File already exists no need to filter")
        filter_job = None

    # grab all relevant genes across all chromosomes
    # simpler if gene details are condensed to one file
    gene_dict: dict[str, dict] = {}
    if chromosomes == "all":
        # autosomes only for now
        chromosomes_list = list(np.arange(22) + 1)
    else:
        chromosomes_list = chromosomes.split(" ")
    for chromosome in chromosomes_list:
        geneloc_tsv_path = dataset_path(
            os.path.join(
                expression_files_prefix,
                "gene_location_files",
                f"GRCh38_geneloc_chr{chromosome}.tsv",
            )
        )

        # concatenating across chromosomes to have a single dict
        gene_dict.update(make_gene_loc_dict(geneloc_tsv_path))

    # isolate to the genes we are interested in
    if genes is not None:
        genes_of_interest = genes.split(" ")
    else:
        genes_of_interest = list(gene_dict.keys())

    # processing cell types (needs to be passed as a single script for click to like it)
    celltype_list = celltypes.split(" ")
    logging.info(f"Cell types to run: {celltype_list}")

    # only run this for relevant genes
    _plink_genes = set()
    for celltype in celltype_list:
        expression_tsv_path = dataset_path(
            os.path.join(
                expression_files_prefix,
                "expression_files",
                f"{celltype}_expression.tsv",
            )
        )
        logging.info(f"before extracting {celltype}-expressed genes - plink files")
        _plink_genes |= set(extract_genes(genes_of_interest, expression_tsv_path))
        logging.info(f"after extracting {celltype}-expressed genes - plink files")
    plink_genes = list(sorted(_plink_genes))
    logging.info(f"Done selecting genes, total number: {len(plink_genes)}")

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.job.Job] = []

    def manage_concurrency_for_job(job: hb.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_gene_concurrency:
            job.depends_on(_dependent_jobs[-max_gene_concurrency])
        _dependent_jobs.append(job)

    # for each gene, extract relevant variants (in window + with some annotation)
    # submit a job for each gene (export genotypes to plink)
    dependencies_dict: Dict[str, hb.job.Job] = {}
    plink_root = output_path("plink_files")
    logging.info("before glob (bim files)")
    storage_client = storage.Client()
    bucket = get_config()["storage"]["default"]["default"].removeprefix("gs://")
    prefix = os.path.join(get_config()["workflow"]["output_prefix"], "plink_files/")

    bim_files = set(
        f"gs://{bucket}/{filepath.name}"
        for filepath in storage_client.list_blobs(bucket, prefix=prefix, delimiter="/")
        if filepath.name.endswith("bim")
    )

    logging.info(f"after glob: {len(bim_files)} bim files already exist")

    for gene in plink_genes:

        # final path for this gene - generate first (check syntax)
        plink_file = os.path.join(plink_root, gene)
        gene_dict[gene]["plink"] = plink_file

        # if the plink output exists, do not re-generate it
        if f"{plink_file}.bim" in bim_files:
            continue

        plink_job = batch.new_python_job(f"Create plink files for: {gene}")
        manage_concurrency_for_job(plink_job)
        copy_common_env(plink_job)
        if filter_job:
            plink_job.depends_on(filter_job)

        plink_job.image(CELLREGMAP_IMAGE)
        plink_job.call(
            get_promoter_variants,
            mt_path=output_mt_path,
            ht_path=anno_ht_path,
            gene_details=gene_dict[gene],
            window_size=window_size,
            plink_file=plink_file,
        )
        dependencies_dict[gene] = plink_job

    # the next phase will be done for each cell type
    for celltype in celltype_list:
        expression_tsv_path = dataset_path(
            os.path.join(
                expression_files_prefix,
                "expression_files",
                f"{celltype}_expression.tsv",
            )
        )
        logging.info(
            f"before extracting {celltype}-expressed genes to run association for"
        )
        genes_list = extract_genes(genes_of_interest, expression_tsv_path)
        logging.info(
            f"after extracting {celltype}-expressed genes: run association for {len(genes_list)} genes"
        )
        # logging.info(f"Genes to run: {genes_list}")
        if not genes_list:
            logging.info("No genes to run, exit!")
            continue

        gene_run_jobs = []

        # cell_type_root = output_path(celltype)
        logging.info(f"before glob: pv files for {celltype}")
        storage_client = storage.Client()
        bucket = get_config()["storage"]["default"]["default"].removeprefix("gs://")
        prefix = f"{get_config()['workflow']['output_prefix']}/{celltype}/"
        existing_files = set(
            f"gs://{bucket}/{filepath.name}"
            for filepath in storage_client.list_blobs(
                bucket, prefix=prefix, delimiter="/"
            )
            if filepath.name.endswith("_pvals.csv")
        )
        # existing_files = list(to_path(cell_type_root).glob("*_pvals.csv"))
        logging.info(f"after glob: {len(existing_files)} pv files for {celltype}")

        for gene in genes_list:

            # wrapped this with output_path
            pv_file = output_path(f"{celltype}/{gene}_pvals.csv")

            # check if running is required
            if pv_file in existing_files:
                logging.info(f"We already ran associations for {gene}!")
                continue

            if gene_dict[gene]["plink"] is None:
                logging.info(f"No plink files for {gene}, exit!")
                continue

            plink_output_prefix = gene_dict[gene]["plink"]
            # prepare input files
            run_job = batch.new_python_job(f"Run association for: {celltype}, {gene}")
            manage_concurrency_for_job(run_job)
            copy_common_env(run_job)
            if dependency := dependencies_dict.get(gene):
                run_job.depends_on(dependency)
            run_job.image(SAIGE_QTL_IMAGE)
            j = batch.new_job("Fit null model")
            j.image("some/R:image")
            step1_cmd = build_fit_null_command()
            j.command(step1_cmd)
            # the python_job.call only returns one object
            # the object is a file containing y_df, geno_df, kinship_df
            # all pickled into a file
            input_results = run_job.call(
                prepare_input_files,
                gene_name=gene,
                cell_type=celltype,
                genotype_file_bed=plink_output_prefix + ".bed",
                genotype_file_bim=plink_output_prefix + ".bim",
                genotype_file_fam=plink_output_prefix + ".fam",
                phenotype_file=expression_tsv_path,
                kinship_file=None,  # change this to work when this file is needed
                sample_mapping_file=sample_mapping_file_tsv,
            )
            # run association in the same python env
            gene_run_jobs.append(run_job)
            run_job.call(
                run_gene_association,
                gene_name=gene,
                prepared_inputs=input_results,
                pv_path=pv_file,
            )

        # combine all p-values across all chromosomes, genes (per cell type)
        summarise_job = batch.new_python_job(f"Summarise all results for {celltype}")
        copy_common_env(summarise_job)
        summarise_job.depends_on(*gene_run_jobs)
        summarise_job.image(MULTIPY_IMAGE)
        pv_all_filename_csv = str(
            output_path(f"{celltype}_all_pvalues.csv", "analysis")
        )
        summarise_job.call(
            summarise_association_results,
            celltype=celltype,
            pv_all_filename_str=str(pv_all_filename_csv),
        )

    # set jobs running
    batch.run(wait=False)


if __name__ == "__main__":
    saige_pipeline()
