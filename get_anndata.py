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
import logging
import math
import hail as hl
import hailtop.batch.job as hb_job
import pandas as pd
import scanpy as sc
from os.path import join
from pathlib import Path
from typing import List

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, init_batch, output_path


def filter_lowly_expressed_genes(expression_adata, min_pct=1) -> sc.AnnData:
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
    gene_info_df_path: str,
    gene: str,
    window_size: int,
    out_path: str,
    chrom_len: int,
):
    """
    Get gene cis window file
    Args:
        gene_info_df_path (str): path to whole adata object
        gene (str): gene name
        window_size (int): bp to consider in window, up and downstream of gene
        out_path (str): path we're writing to
        chrom_len (int): length of chromosome
    """

    gene_info_df = copy_h5ad_local_and_open(gene_info_df_path).var

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
    expression_adata_path: str,
    # sample_covs_df: pd.DataFrame,
    # celltype_covs_df: pd.DataFrame,
    sample_covs_file: str,
    celltype_covs_file: str,
    out_path: str,
    fill_in_sex: bool = True,
    fill_in_age: bool = True,
):
    """
    Combine expression and covariates into a single file
     Args:
        gene (str): gene name
        expression_adata_path (str): path to expression object all genes
        sample_covs_df (pd.DataFrame): sex, age
        celltype_covs_df (pd.DataFrame): celltype specific covs
        out_path (str): path we're writing to
    """
    expression_adata = copy_h5ad_local_and_open(expression_adata_path)

    # open dataframes
    sample_covs_df = pd.read_csv(sample_covs_file)
    sample_covs_df['individual'] = sample_covs_df['sample_id']
    logging.info('sample covariate file opened')
    celltype_covs_df = pd.read_csv(celltype_covs_file, index_col=0)
    logging.info(f'cell covariate file opened')

    # determine average age to fill in later
    if fill_in_age:
        mean_age = sample_covs_df['age'].mean()
    cell_ind_df = expression_adata.obs.loc[:, ['cell', 'individual']]
    sample_covs_cells_df = cell_ind_df.merge(
        sample_covs_df, on='individual', how='left'
    )
    sample_covs_cells_df.index = sample_covs_cells_df['cell']
    # fill missing values for sex and age
    if fill_in_sex:
        sample_covs_cells_df['sex'] = sample_covs_cells_df['sex'].fillna(0)
    if fill_in_age:
        sample_covs_cells_df['age'] = sample_covs_cells_df['age'].fillna(mean_age)
    gene_adata = expression_adata[:, expression_adata.var['gene_name'] == gene]
    expr_df = pd.DataFrame(
        data=gene_adata.X.todense(), index=gene_adata.obs.index, columns=[gene]
    )
    pheno_cov_df = pd.concat([sample_covs_cells_df, expr_df, celltype_covs_df], axis=1)
    with to_path(out_path).open('w') as pcf:
        pheno_cov_df.to_csv(pcf, index=False)


def copy_h5ad_local_and_open(adata_path: str) -> sc.AnnData:
    """
    Copy h5ad file to local, then open it
    Args:
        adata_path (str): path to h5ad file

    Returns:
        sc.AnnData: opened AnnData object
    """
    # copy in h5ad file to local, then load it
    expression_h5ad_path = to_path(adata_path).copy('here.h5ad')
    expression_adata = sc.read_h5ad(expression_h5ad_path)
    assert isinstance(expression_adata, sc.AnnData), type(expression_adata)
    return expression_adata


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--anndata-files-prefix', default='saige-qtl/anndata_objects_from_HPC')
@click.option(
    '--celltype-covs-files-prefix', default='saige-qtl/celltype_covs_from_HPC'
)
@click.option('--sample-covs-files-prefix', default='saige-qtl/input_files/covariates')
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
    # set this up with the default (scanpy) python image
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

    init_batch()

    # sample level covariates (age + sex)
    # age from metamist, sex from somalier + Vlad's file for now
    sample_covs_file = dataset_path(
        f'{sample_covs_files_prefix}/sex_age_tob_bioheart.csv'
    )
    # sample_covs_df = pd.read_csv(sample_covs_file)
    # sample_covs_df['individual'] = sample_covs_df['sample_id']
    # logging.info('sample covariate file opened')

    for celltype in celltypes.split(','):

        # extract cell-level covariates
        # expression PCs, cell type specific
        celltype_covs_file = dataset_path(
            f'{celltype_covs_files_prefix}/{celltype}_expression_pcs.csv'
        )
        # celltype_covs_df = pd.read_csv(celltype_covs_file, index_col=0)
        # logging.info(f'cell covariate for {celltype} file opened')

        for chromosome in chromosomes.split(','):
            chrom_len = hl.get_reference('GRCh38').lengths[chromosome]

            # copy in h5ad file to local, then load it
            expression_h5ad_path = dataset_path(
                f'{anndata_files_prefix}/{celltype}_{chromosome}.h5ad'
            )
            expression_adata = copy_h5ad_local_and_open(expression_h5ad_path)
            logging.info(
                f'AnnData for {celltype}, {chromosome} opened: {expression_adata.shape[1]} total genes'
            )

            # extract genes expressed in at least X% cells
            expression_adata = filter_lowly_expressed_genes(
                expression_adata, min_pct=min_pct_expr
            )
            logging.info(
                f'AnnData for {celltype}, {chromosome} filtered: {expression_adata.shape[1]} genes left'
            )

            # dump the adata to a local file
            tmp_adata_name = f'{celltype}_{chromosome}.h5ad'
            expression_adata.write_h5ad(filename=Path(tmp_adata_name))

            # then write that to GCP
            tmp_path = join(get_config()['storage']['default']['tmp'], tmp_adata_name)
            with open(tmp_adata_name, 'rb') as reader:
                with hl.hadoop_open(tmp_path, 'wb') as writer:
                    for line in reader:
                        writer.write(line)

            # start up some jobs for each gene
            for gene in expression_adata.var['gene_name'][0:10]:

                # make pheno cov file
                pheno_cov_filename = to_path(
                    output_path(
                        f'pheno_cov_files/{celltype}/{chromosome}/{gene}_{celltype}_pheno_cov.csv'
                    )
                )
                if not pheno_cov_filename.exists():
                    pheno_cov_job = get_batch().new_python_job(
                        name=f'pheno cov file: {gene}, {celltype}'
                    )
                    pheno_cov_job.storage('10G')
                    pheno_cov_job.call(
                        make_pheno_cov,
                        gene,
                        str(tmp_path),
                        str(sample_covs_file),
                        str(celltype_covs_file),
                        # sample_covs_df,
                        # celltype_covs_df,
                        str(pheno_cov_filename),
                    )
                    manage_concurrency(pheno_cov_job)
                    logging.info(f'pheno cov job for {gene}, {celltype} scheduled')

                # make cis window file
                gene_cis_filename = to_path(
                    output_path(
                        f'cis_window_files/{chromosome}/{gene}_{cis_window_size}bp.csv'
                    )
                )
                if not gene_cis_filename.exists():
                    gene_cis_job = get_batch().new_python_job(
                        name=f'gene cis file: {gene}'
                    )
                    gene_cis_job.storage('10G')
                    gene_cis_job.call(
                        get_gene_cis_info,
                        str(tmp_path),
                        gene,
                        cis_window_size,
                        str(gene_cis_filename),
                        chrom_len,
                    )
                    manage_concurrency(gene_cis_job)
                    logging.info(f'cis window job for {gene} scheduled')

            ## these aren't necessary with fewer loops, for removal?

            # del expression_adata
            #
            # # delete the local here.h5ad file
            # to_path(expression_h5ad_path).unlink()

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
