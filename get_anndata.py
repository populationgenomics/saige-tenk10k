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


import click
import math
import hail as hl
import hailtop.batch.job as hb_job
import pandas as pd
import scanpy as sc
from typing import List

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, output_path


def filter_lowly_expressed_genes(expression_adata, min_pct=5) -> sc.AnnData:
    """
    Remove genes with low expression across cells

    Args:
        expression_adata (sc.AnnData): adata with all genes
        min_pct (int): min percentage to retain

    Returns:
        filtered AnnData object
    """
    n_all_cells = len(expression_adata.obs.index)
    min_cells = math.ceil((n_all_cells * min_pct) / 100)
    sc.pp.filter_genes(expression_adata, min_cells=min_cells)
    assert isinstance(expression_adata, sc.AnnData), type(expression_adata)

    return expression_adata


def get_gene_cis_info(
    gene_info_df: pd.DataFrame,
    gene: str,
    window_size: int,
    out_path: str,
    chrom_len: int,
):
    """
    Get gene cis window file
    Args:
        gene_info_df (pd.DataFrame): gene chrom, start, end
        gene (str): gene name
        window_size (int): bp to consider in window, up and downstream of gene
        out_path (str): path we're writing to
        chrom_len (int): length of chromosome
    """
    # select the gene from df
    gene_info_gene = gene_info_df[gene_info_df['gene_name'] == gene]
    # get gene chromosome
    chrom = gene_info_gene['chr'][0]
    # get gene body position (start and end) and add window
    left_boundary = max(1, int(gene_info_gene['start'][0]) - window_size)
    right_boundary = min(int(gene_info_gene['end'][0]) + window_size, chrom_len)
    data = {'chromosome': chrom, 'start': left_boundary, 'end': right_boundary}
    gene_cis_df = pd.DataFrame(data, index=[gene])
    with to_path(out_path).open('w') as gcf:
        gene_cis_df.to_csv(gcf, index=False)


def make_pheno_cov(
    gene: str,
    expression_adata: sc.AnnData,
    sample_covs_df: pd.DataFrame,
    celltype_covs_df: pd.DataFrame,
    out_path: str,
):
    """
    Combine expression and covariates into a single file
     Args:
        gene (str): gene name
        expression_adata (sc.AnnData): expression object all genes
        sample_covs_df (pd.DataFrame): sex, age
        celltype_covs_df (pd.DataFrame): celltype specific covs
        out_path (str): path we're writing to
    """
    cell_ind_df = expression_adata.obs['Individual', 'Cell']
    expr_df = expression_adata.X[gene]
    pheno_cov_df = pd.concat(cell_ind_df, expr_df, sample_covs_df, celltype_covs_df)
    with to_path(out_path).open('w') as pcf:
        pheno_cov_df.to_csv(pcf, index=False)


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--anndata-files-prefix', default='saige-qtl/anndata_objects_from_HPC')
@click.option(
    '--celltype-covs-files-prefix', default='saige-qtl/celltype_covs_from_HPC'
)
@click.option('--sample-covs-files-prefix', default='saige-qtl/input_files/covariates/')
@click.option('--min-pct-expr', type=int, default=5)
@click.option('--cis-window-size', type=int, default=100000)
@click.option(
    '--concurrent-job-cap',
    type=int,
    default=100,
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
    celltype_covs_files_prefix: str,
    sample_covs_files_prefix: str,
    min_pct_expr: int,
    cis_window_size: int,
    concurrent_job_cap: int,
):
    """
    Run expression processing pipeline
    """
    # set this up with the default python image
    get_batch(
        default_python_image=get_config()['images']['scanpy'],
        name='prepare all gene files',
    )
    all_jobs: List[hb_job.Job] = []

    def manage_concurrency(new_job: hb_job.Job):
        """
        Manage concurrency, so that there is a cap on simultaneous jobs
        Args:
            new_job (hb_job.Job): a new job to add to the stack
        """
        if len(all_jobs) > concurrent_job_cap:
            new_job.depends_on(all_jobs[-concurrent_job_cap])
        all_jobs.append(new_job)

    hl.init()
    hl.default_reference(hl.get_reference('GRCh38'))

    # extract samples we actually want to test

    # sample level covariates (age + sex)
    # age from metamist, sex from somalier + Vlad's file for now
    sample_covs_file = dataset_path(
        f'{sample_covs_files_prefix}sex_age_tob_bioheart.csv'
    )
    sample_covs_df = pd.read_csv(sample_covs_file)

    for celltype in celltypes.split(','):

        # extract cell-level covariates
        # expression PCs, cell type specific
        celltype_covs_file = dataset_path(
            f'{celltype_covs_files_prefix}/{celltype}_expression_pcs.csv'
        )
        celltype_covs_df = pd.read_csv(celltype_covs_file)

        for chromosome in chromosomes.split(','):
            chrom_len = hl.get_reference('GRCh38').lengths[chromosome]

            # copy in h5ad file to local, then load it
            expression_h5ad_path = to_path(
                dataset_path(f'{anndata_files_prefix}/{celltype}_{chromosome}.h5ad')
            ).copy('here.h5ad')
            expression_adata = sc.read(expression_h5ad_path)
            assert isinstance(expression_adata, sc.AnnData), type(expression_adata)

            # extract genes expressed in at least X% cells
            expression_adata = filter_lowly_expressed_genes(
                expression_adata, min_pct=min_pct_expr
            )

            # start up some jobs for all each gene
            for gene in expression_adata.var['gene_name']:

                # make pheno cov file
                pheno_cov_filename = to_path(
                    output_path(f'expression_files/{gene}_pheno_cov.csv')
                )
                if not pheno_cov_filename.exists():
                    pheno_cov_job = get_batch().new_python_job(name='pheno cov file')
                    pheno_cov_job.call(
                        make_pheno_cov,
                        gene,
                        expression_adata,
                        sample_covs_df,
                        celltype_covs_df,
                    )
                    manage_concurrency(pheno_cov_job)

                # make cis window file
                gene_cis_filename = to_path(
                    output_path(f'cis_window_files/{gene}_{cis_window_size}bp.csv')
                )
                if not gene_cis_filename.exists():
                    gene_cis_job = get_batch().new_python_job(name=f'gene cis file: {gene}')
                    gene_cis_job.call(
                        get_gene_cis_info,
                        expression_adata.var,
                        gene,
                        cis_window_size,
                        str(gene_cis_filename),
                        chrom_len,
                    )
                    manage_concurrency(gene_cis_job)
            del expression_adata

            # delete the local here.h5ad file
            expression_h5ad_path.unlink()

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
