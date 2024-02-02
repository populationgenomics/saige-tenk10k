#!/usr/bin/env python3

"""
This script will

- open anndata expression files
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
    python3 get_anndata.py --celltypes CD4_Naive --chromosomes chr1


"""

from cpg_utils import to_path
from cpg_utils.hail_batch import (
    dataset_path,
    get_batch,
    get_config,
    # init_batch,
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
    expression_adata = sc.pp.filter_genes(expression_adata, min_cells=min_cells)
    assert isinstance(expression_adata, sc.AnnData)

    return expression_adata


# def get_celltype_covariates(
#     expression_files_prefix: str,
#     cell_type: str,
# ):
#     """Obtain cell type specific covariates

#     Input:
#     - cell type of interest
#     - covariate files prefix

#     Output: covariate df for cell type of interest
#     """
#     covs_tsv_path = dataset_path(
#         f'{expression_files_prefix}/expression_pcs/{cell_type}.csv'
#     )
#     covs_df = pd.read_csv(covs_tsv_path, sep=',', index_col=0)
#     return covs_df


def get_gene_cis_info(gene_info_df, gene: str, window_size: int):
    """Get gene cis window file"""
    # select the gene from df
    gene_info_gene = gene_info_df[gene_info_df['gene'] == gene]
    # get gene chromosome
    chrom = gene_info_gene['chr']
    # get gene body position (start and end) and add window
    left_boundary = max(1, int(gene_info_gene['start']) - window_size)
    right_boundary = min(
        int(gene_info_gene['end']) + window_size,
        hl.get_reference('GRCh38').lengths[chrom],
    )
    data = {'chromosome': chrom, 'start': left_boundary, 'end': right_boundary}
    # check if I need an index at all
    return pd.DataFrame(data, index=[gene])


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option(
    '--anndata-files-prefix', default=dataset_path('saige-qtl/anndata_objects_from_HPC')
)
@click.option('--min-pct-expr', type=int, default=5)
@click.option('--cis-window-size', type=int, default=100000)
@click.option(
    '--max-gene-concurrency',
    type=int,
    default=50,
    help=(
        'To avoid resource starvation, set this concurrency to limit '
        'horizontal scale. Higher numbers have a better walltime, but '
        'risk jobs that are stuck (which are expensive)'
    ),
)
def main(
    celltypes: str,
    chromosomes: str,
    anndata_files_prefix: str,
    min_pct_expr: int,
    cis_window_size: int,
    # max_gene_concurrency=int,
):
    """
    Run expression processing pipeline
    """
    # extract samples we actually want to test

    # extract sample level covariates
    # age from metamist
    # sex from somalier

    for celltype in celltypes.split(','):

        # extract cell-level covariates
        # expression PCs, cell type specific

        for chromosome in chromosomes.split(','):
            expression_h5ad_path = to_path(
                dataset_path(f'{anndata_files_prefix}/{celltype}_{chromosome}.h5ad')
            ).copy('here.h5ad')
            expression_adata = sc.read(expression_h5ad_path)

            # extract genes expressed in at least X% cells
            expression_adata = filter_lowly_expressed_genes(
                expression_adata, min_pct=min_pct_expr
            )

            # for each gene
            genes = expression_adata.var['gene_name']

            for gene in genes:
                # get expression
                # make pheno cov file
                # pheno_cov_filename = to_path(
                #     output_path(f'expression_files/{gene}_pheno_cov.csv')
                # )

                # make cis window file
                gene_cis_filename = to_path(
                    output_path(f'cis_window_files/{gene}_{cis_window_size}bp.csv')
                )
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
