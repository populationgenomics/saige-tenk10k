#!/usr/bin/env python3

"""
This script will

- open anndata expression files
- open covariate files
- create pheno_cov_files
- create cis window files

these files will be used as inputs for the
SAIGE-QTL association pipeline.

To run:

analysis-runner \
    --description "make expression input files" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
    --memory "5Gi" --storage "5Gi" \
    python3 get_anndata.py --celltypes CD4_Naive --chromosomes chr1


"""

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    get_batch,
    get_config,
    init_batch,
    output_path,
)
import click
import math
import hail as hl
import pandas as pd
import scanpy as sc

SCANPY_IMAGE = get_config()['images']['scanpy']


def filter_lowly_expressed_genes(expression_adata, min_pct=5) -> sc.AnnData:
    """Remove genes with low expression across cells

    Input: adata with all genes

    Output: adata filtered
    """
    n_all_cells = len(expression_adata.obs.index)
    min_cells = math.ceil((n_all_cells * min_pct) / 100)
    sc.pp.filter_genes(expression_adata, min_cells=min_cells)
    assert isinstance(expression_adata, sc.AnnData), type(expression_adata)

    return expression_adata


def get_gene_cis_info(gene_info_df, gene: str, window_size: int):
    """Get gene cis window file"""
    print(type(gene_info_df))
    print(gene_info_df.head())
    # select the gene from df
    gene_info_gene = gene_info_df[gene_info_df['gene_name'] == gene]
    print(gene_info_gene.head())
    # get gene chromosome
    chrom = gene_info_gene['chr'][0]
    print(type(chrom))
    print(chrom)
    # get gene body position (start and end) and add window
    left_boundary = max(1, int(gene_info_gene['start'][0]) - window_size)
    right_boundary = min(
        int(gene_info_gene['end'][0]) + window_size,
        hl.get_reference('GRCh38').lengths[chrom],
    )
    data = {'chromosome': chrom, 'start': left_boundary, 'end': right_boundary}
    # check if I need an index at all
    return pd.DataFrame(data, index=[gene])


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--anndata-files-prefix', default='saige-qtl/anndata_objects_from_HPC')
# @click.option(
#     '--celltype-covs-files-prefix', default='saige-qtl/celltype_covs_from_HPC'
# )
# @click.option(
#     '--sample-covs-files-prefix', default='saige-qtl/input_files/covariates/'
# )
@click.option('--min-pct-expr', type=int, default=5)
@click.option('--cis-window-size', type=int, default=100000)
# @click.option(
#     '--max-gene-concurrency',
#     type=int,
#     default=50,
#     help=(
#         'To avoid resource starvation, set this concurrency to limit '
#         'horizontal scale. Higher numbers have a better walltime, but '
#         'risk jobs that are stuck (which are expensive)'
#     ),
# )
def main(
    celltypes: str,
    chromosomes: str,
    anndata_files_prefix: str,
    # celltype_covs_files_prefix: str,
    # sample_covs_files_prefix: str,
    min_pct_expr: int,
    cis_window_size: int,
    # max_gene_concurrency=int,
):
    """
    Run expression processing pipeline
    """
    init_batch()
    # batch = get_batch('gene expression processing pipeline')
    # extract samples we actually want to test

    # extract sample level covariates
    # age from metamist
    # sex from somalier
    # sample_covs_file = dataset_path(f'{sample_covs_files_prefix}sex_tob_bioheart.csv')

    for celltype in celltypes.split(','):

        # extract cell-level covariates
        # # expression PCs, cell type specific
        # celltype_covs_file = dataset_path(
        #     f'{celltype_covs_files_prefix}/{celltype}_expression_pcs.csv'
        # )
        # celltype_covs_df = pd.read_csv(celltype_covs_file)

        for chromosome in chromosomes.split(','):
            expression_h5ad_path = to_path(
                dataset_path(f'{anndata_files_prefix}/{celltype}_{chromosome}.h5ad')
            ).copy('here.h5ad')
            expression_adata = sc.read(expression_h5ad_path)
            assert isinstance(expression_adata, sc.AnnData), type(expression_adata)

            # extract genes expressed in at least X% cells
            expression_adata = filter_lowly_expressed_genes(
                expression_adata, min_pct=min_pct_expr
            )
            print(expression_adata.shape)

            # for each gene
            genes = expression_adata.var['gene_name']
            # to test if memory error is due to too many genes, reduce
            genes = genes[0:10]
            # print(genes)

            for gene in genes:
                # print(gene)
                # get expression
                # make pheno cov file
                # pheno_cov_filename = to_path(
                #     output_path(f'expression_files/{gene}_pheno_cov.csv')
                # )

                # make cis window file
                gene_cis_filename = to_path(
                    output_path(f'cis_window_files/{gene}_{cis_window_size}bp.csv')
                )
                # gene_cis_job = batch.new_python_job(name='gene cis file')
                gene_cis_df = get_gene_cis_info(
                    gene_info_df=expression_adata.var,
                    gene=gene,
                    window_size=cis_window_size,
                )
                # write
                with gene_cis_filename.open('w') as gcf:
                    gene_cis_df.to_csv(gcf, index=False)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
