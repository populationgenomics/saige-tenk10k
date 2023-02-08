import os

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    output_path,
)

import pandas as pd
import scanpy as sc
import xarray as xr

# when running in analysis runner, use
# --output-dir "tob_wgs_rv/saige_qtl/"

# region INPUT_FILES

# geneloc files (tsv)
# these are split by chromosome and contain info about the genes
# chrom, ensembl_gene_id, gene_name, strand, start, end
geneloc_files_prefix = 'scrna-seq/grch38_association_files'
chromosome = 1
geneloc_tsv_path = dataset_path(
    os.path.join(
        geneloc_files_prefix,
        "gene_location_files",
        f"GRCh38_geneloc_chr{chromosome}.tsv",
    )
)

# gene expression files (h5ad)
# also one per chromosome, consider whether these are getting too big
# and should be split by cell type also?
expression_files_prefix = 'tob_wgs_rv/saige_qtl/input'
expression_h5ad_path = dataset_path(
    os.path.join(
        expression_files_prefix,
        "expression_objects",
        f"sce{chromosome}.h5ad",
    )
)

# endregion INPUT_FILES

# region FUNCTIONS

def prepare_pheno_cov_file(
    gene_name: str,
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
    pheno_cov file path? maybe no need to return anything, just write to file
    """

    pheno_cov_filename = to_path(
        output_path(f"input/pheno_cov_files/{gene_name}_{cell_type}.tsv")
    )

    # this file will map cells to donors and onek1k ids to cpd ones
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
    phenotype = phenotype.sel(cell=sample_mapping['cell_barcode'].values)

    # delete large files to free up memory
    del mat
    del mat_df


    # read in covariate file (tsv)
    # this file is defined at cell level, as:
    # cell barcode | cov1 | cov2 | ... | cov N
    covs = pd.read_csv(cov_file, sep="\t", index_col=0)

    # this file will map different IDs (and OneK1K ID to CPG ID) as well as donors to cells
    #         CPG ID  | OneK1K ID |    cell barcode
    # e.g.,   CPG7385 |  686_687  | AAACCTGCAACGATCT-1
    sample_mapping = pd.read_csv(dataset_path(sample_mapping_file), sep="\t")

    # add individual ID to covariates
    cov_samples = covs.merge(sample_mapping, on='cell_barcode')


    # phenotype
    # select gene
    y = phenotype.sel(gene=gene_name)
    # y = quantile_gaussianize(y)  # do not do this here (will use a Poisson likelihood)
    del phenotype  # delete to free up memory
    # make data frame to save as tsv
    y_df = pd.DataFrame(
        data=y.values.reshape(y.shape[0], 1), index=y.sample.values, columns=[gene_name]
    )

    # make final data frame
    # columns = y | cov1 | ... | covN | indID
    pheno_cov_df = y_df.merge(cov_samples, on="barcode")
    print(pheno_cov_df.head())

    # save files
    with pheno_cov_filename.open("w") as pcf:
        pheno_cov_df.to_csv(pcf, index=False, sep="\t")


# endregion FUNCTIONS