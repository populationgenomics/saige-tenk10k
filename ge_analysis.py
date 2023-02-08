import os
import logging

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

# region MISC

# make gene loc info into a dict with genes as keys
def make_gene_loc_dict(file) -> dict[str, dict]:
    """
    Turn gene information into a dictionary
    to avoid opening the gene loc file for every gene
    """
    from csv import DictReader

    gene_dict = {}

    with open(to_path(file)) as handle:
        reader = DictReader(handle, delimiter="\t")

        for row in reader:
            gene_dict[row["gene_name"]] = row

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

    logging.info(f"Total genes to run: {len(list(sorted(genes)))}")

    return list(sorted(genes))

# endregion MISC

# region PHENO_COV_FILE

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
    covs = pd.read_csv(cov_file, sep="\t", index_col=0)

    # add individual ID to covariates (this will also subset covs to right cells)
    cov_samples = covs.merge(sample_mapping, on='cell_barcode')

    # phenotype
    # select gene
    expr = phenotype.sel(gene=gene_name)
    del phenotype  # delete to free up memory
    # make data frame to save as tsv
    expr_df = pd.DataFrame(
        data=expr.values.reshape(expr.shape[0], 1), index=expr.sample.values, columns=[gene_name]
    )

    # make final data frame
    # columns = y | cov1 | ... | covN | indID
    pheno_cov_df = expr_df.merge(cov_samples, on="cell_barcode")

    # save file
    with pheno_cov_filename.open("w") as pcf:
        pheno_cov_df.to_csv(pcf, index=False, sep="\t")


# endregion PHENO_COV_FILE