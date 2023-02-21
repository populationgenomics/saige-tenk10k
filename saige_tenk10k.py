#!/usr/bin/env python3
# pylint: disable=import-error,import-outside-toplevel,missing-function-docstring,missing-module-docstring,no-value-for-parameter,too-many-arguments,too-many-locals,wrong-import-order,wrong-import-position,consider-using-enumerate,chained-comparison

__author__ = 'annacuomo'

"""
Hail Batch workflow for the rare-variant association analysis, including:

- get relevant variants around a gene and export genotypes as plink files,
- generate other input files for association tests (phenotype, covariates, groups),
- run association tests.

More details in README
output files in tob_wgs_genetics/saige_qtl/output
"""

# import python modules
import os
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
import scanpy as sc
import xarray as xr

import hail as hl
import hailtop.batch as hb


# use logging to print statements, display at info level
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr,
)

DEFAULT_JOINT_CALL_MT = dataset_path('mt/v7.mt')
DEFAULT_ANNOTATION_HT = dataset_path(
    'tob_wgs_vep/104/vep104.3_GRCh38.ht'
)  # atm VEP only - add open chromatin info

HAIL_IMAGE = get_config()['workflow']['driver_image']
SAIGE_QTL_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/saige-qtl'
MULTIPY_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16'

# region SUBSET_VARIANTS


# only needs to be run once for a given cohort (e.g., OneK1K)
def filter_variants(
    mt_path: str,  # 'mt/v7.mt'
    samples: list[str],
    output_rv_mt_path: str,  # 'tob_wgs/densified_rv_only.mt'
    output_cv_mt_path: str,  # 'tob_wgs/densified_cv_only.mt'
    vre_plink_path: str,  # 'tob_wgs/vr_plink_2000_variants
    cv_maf_threshold: float = 0.01,
    rv_maf_threshold: float = 0.05,
    vre_mac_threshold: int = 20,
    vre_n_markers: int = 2000,
):
    """Subset hail matrix table

    Input:
    - joint call hail matrix table
    - set of samples for which we have scRNA-seq data
    - file paths for ouputs
    - MAF thresholds to define common / rare variants

    Output 1&2:
    subset hail matrix table, containing only variants that:
    i) are not ref-only, ii) biallelic, iii) meet QC filters,
    and samples that are contained in sc sample list.

    Then, in output1 variants are also rare (freq < rv_maf_threshold)
    and in output that are common (freq > cv_maf_threshold)

    Output 3:
    plink file containing a random subset of 2,000 variants that satisfy i),ii),iii)
    that are additionally sufficiently common (MAC>20) and not in LD
    """
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # subset to relevant samples (samples we have scRNA-seq data for)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # densify
    mt = hl.experimental.densify(mt)

    # filter out low quality variants and consider biallelic SNPs only (no multi-allelic, no ref-only, no indels)
    mt = mt.filter_rows(  # check these filters!
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic (revisit)
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (revisit)
    )

    mt = hl.variant_qc(mt)

    # subset variants for variance ratio estimation
    # minor allele count (MAC) > 20
    tot_counts = mt.variant_qc.AC.sum()  # ????
    vre_mt = mt.filter_rows(
        (mt.variant_qc.AC[1] > vre_mac_threshold) & (mt.variant_qc.AC[1] < tot_counts)
        | (mt.variant_qc.AF[1] < (tot_counts - vre_mac_threshold))
        & (mt.variant_qc.AC[1] > 0)
    )
    # perform LD pruning
    vre_mt = vre_mt.sample_rows(
        p=0.01
    )  # in case this is very costly, subset first a bit
    pruned_variant_table = hl.ld_prune(vre_mt.GT, r2=0.2, bp_window_size=500000)
    vre_mt = vre_mt.filter_rows(hl.is_defined(pruned_variant_table[vre_mt.row_key]))
    # randomly sample {vre_n_markers} variants
    vre_mt = vre_mt.sample_rows((vre_n_markers * 1.1) / vre_mt.count[0])
    vre_mt = vre_mt.head(vre_n_markers)

    # export to plink common variants only for sparse GRM
    from hail.methods import export_plink

    export_plink(vre_mt, vre_plink_path, ind_id=vre_mt.s)

    # filter common variants for single-variant association
    cv_mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > cv_maf_threshold) & (mt.variant_qc.AF[1] < 1)
        | (mt.variant_qc.AF[1] < (1 - cv_maf_threshold)) & (mt.variant_qc.AF[1] > 0)
    )
    cv_mt.write(output_cv_mt_path, overwrite=True)
    logging.info(
        f'Number of common (freq>{cv_maf_threshold*100}%) and QCed biallelic SNPs: {cv_mt.count()[0]}'
    )

    # filter rare variants only (MAF < 5%)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] < rv_maf_threshold) & (mt.variant_qc.AF[1] > 0)
        | (mt.variant_qc.AF[1] > (1 - rv_maf_threshold)) & (mt.variant_qc.AF[1] < 1)
    )
    mt.write(output_rv_mt_path, overwrite=True)
    logging.info(
        f'Number of rare (freq<{rv_maf_threshold*100}%) and QCed biallelic SNPs: {mt.count()[0]}'
    )


# endregion SUBSET_VARIANTS

# region PREPARE_PHENO_COV_FILE


def prepare_pheno_cov_file(
    chromosome: str,
    cell_type: str,
    phenotype_file: str,
    cov_file: str,
    sample_mapping_file: str,
):
    """Prepare pheno+cov file for SAIGE-QTL

    Input:
    phenotype: gene expression (h5ad)
    covariates: cell-level (tsv)

    Output:
    pheno_cov file path (not actually returned, just written to inside function)
    """

    pheno_cov_filename = to_path(
        output_path(f'input/pheno_cov_files/{cell_type}_chr{chromosome}_allgenes.tsv')
    )

    # this file will map different IDs (and OneK1K ID to CPG ID) as well as donors to cells
    #         CPG ID  | OneK1K ID |    cell barcode     | cell type
    # e.g.,   CPG7385 |  686_687  | AAACCTGCAACGATCT-1  |   CD4_T
    sample_mapping = pd.read_csv(
        sample_mapping_file,
        dtype={
            'onek1k_id_long': str,
            'onek1k_id_short': str,
            'cpg_id': str,
            'cell_barcode': str,
            'celltype_label': str,
        },
        index_col=0,
    )
    # subset to relevant cells (given cell_type)
    sample_mapping = sample_mapping['celltype_label' == cell_type]

    # read in phenotype file (scanpy object AnnData)
    # open anndata
    adata = sc.read(phenotype_file)
    # sparse to dense
    mat = adata.raw.X.todense()
    # make pandas dataframe
    mat_df = pd.DataFrame(
        data=mat.T, index=adata.raw.var.index, columns=adata.obs.index
    )
    # turn into xr array
    phenotype = xr.DataArray(
        mat_df.values,
        dims=['trait', 'cell'],
        coords={'trait': mat_df.index.values, 'cell': mat_df.columns.values},
    )
    # consider only correct cells
    phenotype = phenotype.sel(cell=sample_mapping['cell_barcode'].values)

    # delete large files to free up memory
    del mat
    del mat_df

    # read in covariate file (tsv)
    # this file is defined at cell level, as:
    # cell barcode | cov1 | cov2 | ... | cov N
    covs = pd.read_csv(cov_file, sep='\t', index_col=0)

    # add individual ID to covariates (this will also subset covs to right cells)
    cov_samples = covs.merge(sample_mapping, on='cell_barcode')

    # make data frame to save as tsv
    expr_df = pd.DataFrame(
        data=phenotype.values.T,
        index=phenotype.sample.values,
        columns=phenotype.trait.values,
    )

    del phenotype  # delete to free up memory

    # make final data frame
    # columns = y | cov1 | ... | covN | indID
    pheno_cov_df = expr_df.merge(cov_samples, on='cell_barcode')

    # save files
    with pheno_cov_filename.open('w') as pcf:
        pheno_cov_df.to_csv(pcf, index=False, sep='\t')


# endregion PREPARE_PHENO_COV_FILE

# region GET_GENE_SPECIFIC_VARIANTS

def get_promoter_variants(
    mt_path: str,  # output path from function above
    vep_ht_path: str,  
    open_chr_ht_path: str, 
    gene_details: dict[str, str],  # output of make_gene_loc_dict
    window_size: int,  
    plink_file: str,  # 'tob_wgs_rv/saige_qtl/input/plink_files/GENE'
):
    """Subset hail matrix table

    Input:
    mt_path: path to already subsetted hail matrix table
    vep_ht_path: path to VEP HT
    open_chr_ht_path: path to open chromatin HT
    gene_details: dict of info for current gene
    window_size: int, size of flanking region around genes
    plink_file: str, file prefix for writing plink data

    Output:
    Retained variants, that are: 1) regulatory based on annotations
    2) within 50kb up or down-stream of the gene body (or in the gene body itself)
    (on top of all filters done above)

    returns nothing (simply writes out to plink)
    """

    # read hail matrix table object (pre-filtered)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    gene_name = gene_details['gene_name']

    # get relevant chromosome
    chrom = gene_details['chr']

    # subset to genomic window around the gene
    # get gene body position (start and end) and build interval
    left_boundary = max(1, int(gene_details['start']) - window_size)
    right_boundary = min(
        int(gene_details['end']) + window_size,
        hl.get_reference('GRCh38').lengths[chrom],
    )
    # get gene-specific genomic interval
    gene_interval = f'{chrom}:{left_boundary}-{right_boundary}'
    logging.info(f'Interval considered: {gene_interval}')  # 'chr22:23219960-23348287'

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome='GRCh38')]
    )
    mt_path = output_path(f'{gene_name}_in_window.mt', 'tmp')
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f'Number of variants within interval: {mt.count()[0]}')

    # add anotations
    # annotate using VEP
    vep_anno_ht = hl.read_table(vep_ht_path)
    mt = mt.annotate_rows(vep=vep_anno_ht[mt.row_key].vep)
    # annotate using open chromatin info
    oc_anno_ht = hl.read_table(open_chromatin_ht_path)
    mt = mt.annotate_rows(atac=oc_anno_ht[mt.row_key].open_chromatin)

    # filter variants found to be in promoter regions - change?
    mt = mt.filter_rows(
        mt.vep.regulatory_feature_consequences['biotype'].contains('promoter')
    )
    # add aditional filtering for all regions that are open
    # these flags will need to be different for each cell type
    # but can be dealth with upon running (setting annotations in a way that only open regions for that chromatin are run)
    mt = mt.filter_rows(mt.atac.open_chromatin > 0)
    anno_path = output_path(f'{gene_name}_open_promoter_variants.mt', 'tmp')
    mt = mt.checkpoint(anno_path, overwrite=True)  # checkpoint
    logging.info(
        f'Number of rare (freq<5%) QC-passing, biallelic SNPs in promoter and open regions: {mt.count()[0]}'
    )

    # export this as a Hail table for downstream analysis
    # consider exporting the whole MT?
    ht_path = output_path(
        f'summary_hts/{gene_name}_rare_promoter_open_summary.ht', 'analysis'
    )
    ht = mt.rows()
    ht.write(ht_path, overwrite=True)

    # write as group file
    df0 = pd.DataFrame([[gene_name, 'var'],[gene_name, 'anno']])
    variants = ht.locus.collect()
    annos = ht.atac.collect()
    data = [[variants], [annos]]
    df = df0.concat(pd.DataFrame(data), axis=1)
    group_filename = output_path(
        f'input/group_files/{gene_name}_groups.tsv',
    )
    df.to_csv(group_filename)

    # export MT object to PLINK (promoter variants)
    # pylint: disable=import-outside-toplevel
    from hail.methods import export_plink

    export_plink(mt, plink_file, ind_id=mt.s)


# endregion GET_GENE_SPECIFIC_VARIANTS


# region GET_SAIGE_COMMANDS

# Fit null model
def build_fit_null_command(
    pheno_file: str,
    cov_col_list: str,  # PC1
    sample_id_pheno: str,  # IND_ID
    plink_path: str,
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
    saige_command_step1 += f' --useGRMtoFitNULL=FALSE'
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
    plink_prefix: str,
    saige_output_file: str,  # should end in txt? just path?
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


# endregion GET_SAIGE_COMMANDS


# region AGGREGATE_RESULTS


def summarise_association_results(
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

    logging.info('before glob (pv files) - summarise job')
    storage_client = storage.Client()
    bucket = get_config()['storage']['default']['default'].removeprefix('gs://')
    prefix = f"{get_config()['workflow']['output_prefix']}/{celltype}/"
    existing_pv_files = set(
        f'gs://{bucket}/{filepath.name}'
        for filepath in storage_client.list_blobs(bucket, prefix=prefix, delimiter='/')
        if filepath.name.endswith('_results.tsv')
    )
    logging.info(f'after glob - {len(existing_pv_files)} pv files to summarise')

    if len(existing_pv_files) == 0:
        raise Exception('No PV files, nothing to do')

    pv_all_df = pd.concat(
        [
            pd.read_csv(to_path(pv_df), index_col=0, sep='\t')
            for pv_df in existing_pv_files
        ]
    )

    # run qvalues for all tests (multiple testing correction)
    _, qvals = qvalue(pv_all_df['Pvalue'])
    pv_all_df['Qvalue'] = list(qvals)
    _, qvals = qvalue(pv_all_df['Pvalue_Burden'])
    pv_all_df['Qvalue_Burden'] = list(qvals)
    _, qvals = qvalue(pv_all_df['Pvalue_SKAT'])
    pv_all_df['Qvalue_SKAT'] = list(qvals)

    pv_all_filename = to_path(pv_all_filename_str)
    logging.info(f'Write summary results to {pv_all_filename}')
    with pv_all_filename.open('w') as pf:
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
        reader = DictReader(handle, delimiter='\t')

        for row in reader:
            gene_dict[row['gene_name']] = row

    return gene_dict


def extract_genes(gene_list, expression_h5ad_path) -> list[str]:
    """
    Takes a list of all genes and subsets to only those
    present in the expression file of interest
    """
    adata = sc.read(to_path(expression_h5ad_path))
    # consider adding extra filters on expression here
    gene_ids = set(list(adata.raw.var.index))
    genes = set(gene_list).intersection(gene_ids)

    logging.info(f'Total genes to run: {len(list(sorted(genes)))}')

    return list(sorted(genes))


def remove_sc_outliers(df, outliers=None):
    """
    Remove outlier samples, as identified by single-cell analysis
    """
    if outliers is None:
        outliers = ['966_967', '88_88']
    else:
        outliers = outliers.extend(['966_967', '88_88'])
    df = df[-df['OneK1K_ID'].isin(outliers)]

    return df


# endregion MISCELLANEOUS

config = get_config()


@click.command()
@click.option('--celltypes')
@click.option('--geneloc-files-prefix', default='scrna-seq/grch38_association_files')
@click.option('--input-files-prefix', default='tob_wgs_genetics/saige_qtl/input')
@click.option(
    '--sample-mapping-file-tsv',
    default='scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv',  # this needs to be updated
)
@click.option('--mt-path', default=DEFAULT_JOINT_CALL_MT)
@click.option(
    '--anno-ht-path', default=DEFAULT_ANNOTATION_HT
)  # this needs to be updated
@click.option(
    '--chromosomes',
    help='List of chromosome numbers to run rare variant association analysis on. '
    'Space separated, as one argument (Default: all)',
)
@click.option('--genes', default=None)
@click.option('--window_size', default=50000, help='cis size of window around gene')
@click.option(
    '--max-gene-concurrency',
    type=int,
    default=50,
    help=(
        'To avoid resource starvation, set this concurrency to limit horizontal scale. '
        'Higher numbers have a better walltime, but risk jobs that are stuck (which are expensive)'
    ),
)
def saige_pipeline(
    celltypes: str,
    geneloc_files_prefix: str,
    input_files_prefix: str,
    sample_mapping_file_tsv: str,
    mt_path: str,
    anno_ht_path: str,
    chromosomes: str = 'all',
    genes: str | None = None,
    window_size: int = 50000,
    max_gene_concurrency=100,
):

    sb = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch('SAIGE-QTL pipeline', backend=sb)

    # extract individuals for which we have single-cell (sc) data
    sample_mapping_file = pd.read_csv(dataset_path(sample_mapping_file_tsv), sep='\t')
    # we may want to exclude these from the smf directly
    sample_mapping_file = remove_sc_outliers(sample_mapping_file)
    # check column names - CPG_ID would be better?
    sc_samples = sample_mapping_file['InternalID'].unique()

    # filter to QC-passing, rare, biallelic variants
    output_rv_mt_path = output_path('densified_rv_only.mt')
    output_cv_mt_path = output_path('densified_cv_only.mt')
    vre_plink_path = output_path('vr_plink_20k_variants')
    if not to_path(output_rv_mt_path).exists():

        filter_job = batch.new_python_job(name='MT filter job')
        copy_common_env(filter_job)
        filter_job.image(HAIL_IMAGE)
        filter_job.call(
            filter_variants,
            mt_path=mt_path,
            samples=list(sc_samples),
            output_rv_mt_path=output_rv_mt_path,
            output_cv_mt_path=output_cv_mt_path,
            vre_plink_path=vre_plink_path,
        )

    else:
        logging.info('File already exists no need to filter')
        filter_job = None

    # grab all relevant genes across all chromosomes
    # simpler if gene details are condensed to one file
    gene_dict: dict[str, dict] = {}
    if chromosomes == 'all':
        # autosomes only for now
        chromosomes_list = list(np.arange(22) + 1)
    else:
        chromosomes_list = chromosomes.split(' ')
    for chromosome in chromosomes_list:
        # check that these are still appropriate
        geneloc_tsv_path = dataset_path(
            os.path.join(
                geneloc_files_prefix,
                'gene_location_files',
                f'GRCh38_geneloc_chr{chromosome}.tsv',
            )
        )
        # concatenating across chromosomes to have a single dict
        gene_dict.update(make_gene_loc_dict(geneloc_tsv_path))

    # isolate to the genes we are interested in
    if genes is not None:
        genes_of_interest = genes.split(' ')
    else:
        genes_of_interest = list(gene_dict.keys())

    # only run if there is sufficient expression
    _plink_genes = set()
    for chromosome in chromosomes_list:
        expression_h5ad_path = dataset_path(
            os.path.join(
                input_files_prefix,
                'expression_objects',
                f'sce{chromosome}.h5ad',
            )
        )
        logging.info(f'before extracting chrom {chromosome} genes - plink files')
        _plink_genes |= set(extract_genes(genes_of_interest, expression_h5ad_path))
        logging.info(f'after extracting chrom {chromosome} genes - plink files')
    plink_genes = list(sorted(_plink_genes))
    logging.info(f'Done selecting genes, total number: {len(plink_genes)}')

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
    plink_root = output_path('plink_files')
    logging.info('before glob (bim files)')
    storage_client = storage.Client()
    bucket = get_config()['storage']['default']['default'].removeprefix('gs://')
    prefix = os.path.join(get_config()['workflow']['output_prefix'], 'plink_files/')

    bim_files = set(
        f'gs://{bucket}/{filepath.name}'
        for filepath in storage_client.list_blobs(bucket, prefix=prefix, delimiter='/')
        if filepath.name.endswith('bim')
    )

    logging.info(f'after glob: {len(bim_files)} bim files already exist')

    for gene in plink_genes:

        # final path for this gene - generate first (check syntax)
        plink_file = os.path.join(plink_root, gene)
        gene_dict[gene]['plink'] = plink_file

        # if the plink output exists, do not re-generate it
        if f'{plink_file}.bim' in bim_files:
            continue

        plink_job = batch.new_python_job(f'Create plink files for: {gene}')
        manage_concurrency_for_job(plink_job)
        copy_common_env(plink_job)
        if filter_job:
            plink_job.depends_on(filter_job)

        plink_job.image(HAIL_IMAGE)
        plink_job.call(
            get_promoter_variants,
            mt_path=output_rv_mt_path,
            ht_path=anno_ht_path,
            gene_details=gene_dict[gene],
            window_size=window_size,
            plink_file=plink_file,
        )
        dependencies_dict[gene] = plink_job

    # processing cell types (needs to be passed as a single script for click to like it)
    celltype_list = celltypes.split(' ')
    logging.info(f'Cell types to run: {celltype_list}')

    # the next phase will be done for each cell type
    for celltype in celltype_list:
        for chromosome in chromosomes_list:
            # input: gene ID, phenotype file, covariate file
            # output: pheno_cov file
            pheno_cov_job = batch.new_python_job(
                f'Make pheno cov file for: {gene}, {celltype}'
            )
            manage_concurrency_for_job(pheno_cov_job)
            copy_common_env(pheno_cov_job)
            # does not depend on any other jobs
            pheno_cov_job.image(HAIL_IMAGE)
            pheno_cov_job.call(
                prepare_pheno_cov_file,
                gene_name=gene,
                cell_type=celltype,
                phenotype_file=f'sce{chromosome}.h5ad',
                cov_file=f'{celltype}_covs.tsv',
                sample_mapping_file=sample_mapping_file,
            )

    for celltype in celltype_list:
        gene_run_jobs = []
        # need to remind myself of what was happening here
        logging.info(f'before glob: result files for {celltype}')
        storage_client = storage.Client()
        bucket = get_config()['storage']['default']['default'].removeprefix('gs://')
        prefix = f"{get_config()['workflow']['output_prefix']}/{celltype}/"
        existing_files = set(
            f'gs://{bucket}/{filepath.name}'
            for filepath in storage_client.list_blobs(
                bucket, prefix=prefix, delimiter='/'
            )
            if filepath.name.endswith('_results.csv')
        )
        logging.info(f'after glob: {len(existing_files)} result files for {celltype}')
        # end of confusing part
        for gene in plink_genes:

            # wrapped this with output_path
            res_file = output_path(f'output/{celltype}/{gene}_results.csv')

            # check if running is required
            if res_file in existing_files:
                logging.info(f'We already ran associations for {gene}!')
                continue

            if gene_dict[gene]['plink'] is None:
                logging.info(f'No plink files for {gene}, exit!')
                continue

            plink_output_prefix = gene_dict[gene]['plink']

            # input: pheno_cov file, subset plink files
            # output: null model object, variance ratio (VR) estimate file
            fit_null_job = batch.new_job(
                f'Fit null model for: {gene}, {celltype}'
            )
            manage_concurrency_for_job(fit_null_job)
            copy_common_env(fit_null_job)
            # syntax below probably does not work
            dependencies = [filter_job, pheno_cov_job]
            fit_null_job.depends_on(dependencies)
            fit_null_job.image(SAIGE_QTL_IMAGE)
            pheno_cov_filename = output_path(
                f'input/pheno_cov_files/{celltype}_chr{chromosome}_allgenes.tsv'
            )
            # python job creating Rscript command
            cmd = build_fit_null_command(
                pheno_file=pheno_cov_filename,
                cov_col_list='PC1',  # define these when making cov file
                sample_id_pheno='cpg_id',  # check
                plink_path=vre_plink_path,  # this is just for variance ratio estimation?
                output_prefix=output_path('output/saige_model'),
                pheno_col=gene,
                trait_type='count',  # to evaluate vs 'count_NB'
            )
            # regular job submitting the Rscript command to bash
            fit_null_job.command(cmd)

            # input: null model object, VR file, gene specific genotypes (w open chromatin flags)
            # output: summary stats
            run_association_job = batch.new_job(
                f'Run gene set association for: {gene}, {celltype}'
            )
            dependencies = [fit_null_job, plink_job]
            run_association_job.depends_on(dependencies)
            run_association_job.image(SAIGE_QTL_IMAGE)
            # build Rscript command
            cmd = build_run_set_test_command(
                plink_prefix=plink_output_prefix,
                saige_output_file=output_path('output/saige'),
                chrom=chromosome,
                gmmat_model_path=output_path('output/saige_model.rda'),
                variance_ratio_path=output_path('output/saige_model.varianceRatio.txt'),
                group_annotation=f'open_chromatin_{celltype}',
                group_file=dataset_path(f'input/group_files/{gene}_group.txt'),
                allele_order='ref-first',
            )
            # regular job submitting the Rscript command to bash
            run_association_job.command(cmd)

        # combine all p-values across all chromosomes, genes (per cell type)
        summarise_job = batch.new_python_job(f'Summarise all results for {celltype}')
        copy_common_env(summarise_job)
        summarise_job.depends_on(*gene_run_jobs)
        summarise_job.image(MULTIPY_IMAGE)
        pv_all_filename_csv = str(
            output_path(f'{celltype}_all_pvalues.csv', 'analysis')
        )
        summarise_job.call(
            summarise_association_results,
            celltype=celltype,
            pv_all_filename_str=str(pv_all_filename_csv),
        )

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    saige_pipeline()
