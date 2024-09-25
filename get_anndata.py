#!/usr/bin/env python3

"""
This script will

- open anndata expression files
- open cell and samples covariate files
- create pheno_cov_files
- create cis window files

these files will be used as inputs for the
SAIGE-QTL association pipelines (both common and rare).

To run:

analysis-runner \
   --description "make expression input files" \
   --dataset "bioheart" \
   --access-level "full" \
   --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/input_files/240920" \
   --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
   python3 get_anndata.py --celltypes B_naive --chromosomes chr2 \
   --anndata-files-prefix gs://cpg-bioheart-main/saige-qtl/240-libraries/anndata_objects_from_HPC \
   --celltype-covs-files-prefix gs://cpg-bioheart-main/saige-qtl/240-libraries/celltype_covs_from_HPC \
   --sample-covs-file gs://cpg-bioheart-main-analysis/saige-qtl/input_files/240920/covariates/sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv \
   --pc-job-mem=8G

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
import subprocess

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path


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
        out_path (str): path we're writing to (TSV)
        chrom_len (int): length of chromosome
    """

    gene_info_df = copy_h5ad_local_and_open(gene_info_df_path).var

    # select the gene from df
    gene_info_gene = gene_info_df[gene_info_df.index == gene]
    # get gene chromosome
    chrom = gene_info_gene['chr'][0]
    # remove "chr"
    chrom = chrom.replace('chr', '')
    # get gene body position (start and end) and add window
    left_boundary = max(1, int(gene_info_gene['start'][0]) - window_size)
    right_boundary = min(int(gene_info_gene['end'][0]) + window_size, chrom_len)
    data = {'chromosome': chrom, 'start': left_boundary, 'end': right_boundary}
    gene_cis_df = pd.DataFrame(data, index=[gene])
    with to_path(out_path).open('w') as gcf:
        gene_cis_df.to_csv(gcf, index=False, header=False, sep='\t')


def make_pheno_cov(
    gene: str,
    expression_adata_path: str,
    sample_covs_file: str,
    celltype_covs_file: str,
    out_path: str,
):
    """
    Combine expression and covariates into a single file
     Args:
        gene (str): gene name
        expression_adata_path (str): path to expression object all genes
        sample_covs_df (pd.DataFrame): sex, age, genotype PCs
        celltype_covs_df (pd.DataFrame): celltype specific covs
        out_path (str): path we're writing to (TSV)
    """
    expression_adata = copy_h5ad_local_and_open(expression_adata_path)

    # barcoding discrepancy - to be fixed in the next freeze
    expression_adata.obs.index = [
        cell.split("-")[0] for cell in expression_adata.obs.index
    ]
    expression_adata.obs.cell = [
        cell.split("-")[0] for cell in expression_adata.obs.index
    ]

    # open dataframes
    sample_covs_df = pd.read_csv(sample_covs_file)
    sample_covs_df['individual'] = sample_covs_df['sample_id']
    logging.info('sample covariate file opened')
    celltype_covs_df = pd.read_csv(celltype_covs_file, index_col=0)
    logging.info('cell covariate file opened')

    cell_ind_df = expression_adata.obs.loc[
        :, ['cell', 'individual', 'total_counts', 'sequencing_library', 'cohort']
    ]
    # make sequencing_library from categorical to dummy numerical covs
    seq_lib_df = pd.get_dummies(cell_ind_df['sequencing_library']).astype(int)
    # do the same for cohort
    cohort_df = pd.get_dummies(cell_ind_df['cohort']).astype(int)
    cell_ind_df = pd.concat([cell_ind_df, cohort_df, seq_lib_df], axis=1)
    # merge cell and sample covs
    sample_covs_cells_df = cell_ind_df.merge(
        sample_covs_df, on='individual', how='inner'
    )
    sample_covs_cells_df.index = sample_covs_cells_df['cell']
    # drop rows with missing values (SAIGE throws an error otherwise:  https://batch.hail.populationgenomics.org.au/batches/435978/jobs/91)
    sample_covs_cells_df = sample_covs_cells_df.dropna()
    gene_adata = expression_adata[:, expression_adata.var.index == gene]
    gene_name = gene.replace("-", "_")
    expr_df = pd.DataFrame(
        data=gene_adata.X.todense(), index=gene_adata.obs.index, columns=[gene_name]
    )
    # move index (barcode) into a 'cell' column and reset the index - required prior to merging
    # TO DO adjust when final data comes (see issue #97)
    expr_df['cell'] = expr_df.index
    expr_df = expr_df.reset_index(drop=True)
    sample_covs_cells_df = sample_covs_cells_df.reset_index(drop=True)
    celltype_covs_df['cell'] = celltype_covs_df.index
    celltype_covs_df = celltype_covs_df.reset_index(drop=True)

    sample_covs_cells_df = pd.merge(sample_covs_cells_df, expr_df, on='cell')
    pheno_cov_df = pd.merge(
        sample_covs_cells_df, celltype_covs_df, on='cell', how='inner'
    )

    with to_path(out_path).open('w') as pcf:
        pheno_cov_df.to_csv(pcf, index=False, sep='\t')


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
@click.option('--anndata-files-prefix', default='GCP path to anndata file directory')
@click.option(
    '--celltype-covs-files-prefix',
    help='GCP path to celltype covariates file directory',
)
@click.option('--sample-covs-file', help='GCP path to sample covariates file')
@click.option('--min-pct-expr', type=int, default=1)
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
@click.option(
    '--pc-job-cpu',
    type=float,
    default=1,
    help='CPU for each pheno cov job',
)
@click.option(
    '--pc-job-mem',
    type=str,
    default='standard',
    help='Memory for each pheno cov job',
)
@click.option(
    '--cis-job-cpu',
    type=float,
    default=0.25,
    help='CPU for each cis window job',
)
def main(
    celltypes: str,
    chromosomes: str,
    anndata_files_prefix: str,
    celltype_covs_files_prefix: str,
    sample_covs_file: str,
    min_pct_expr: int,
    cis_window_size: int,
    concurrent_job_cap: int,
    pc_job_cpu: float,
    pc_job_mem: str,
    cis_job_cpu: float,
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

    # sample level covariates (age + sex + genotype PCs)
    # age from metamist, sex from somalier + Vlad's file for now

    for celltype in celltypes.split(','):
        # extract cell-level covariates
        # expression PCs, cell type specific
        celltype_covs_file = (
            f'{celltype_covs_files_prefix}/{celltype}_expression_pcs.csv'
        )

        for chromosome in chromosomes.split(','):
            chrom_len = hl.get_reference('GRCh38').lengths[chromosome]

            # copy in h5ad file to local, then load it
            expression_h5ad_path = (
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
            cmd = ["gsutil", "cp", tmp_adata_name, tmp_path]
            subprocess.run(cmd, check=True)

            # start up some jobs for each gene
            for gene in expression_adata.var.index:
                # change hyphens to underscore for R usage
                gene_name = gene.replace("-", "_")
                # make pheno cov file
                pheno_cov_filename = to_path(
                    output_path(
                        f'pheno_cov_files/{celltype}/{chromosome}/{gene_name}_{celltype}_pheno_cov.tsv'
                    )
                )
                if not pheno_cov_filename.exists():
                    pheno_cov_job = get_batch().new_python_job(
                        name=f'pheno cov file: {gene}, {celltype}'
                    )
                    pheno_cov_job.cpu(pc_job_cpu)
                    pheno_cov_job.memory(pc_job_mem)
                    pheno_cov_job.call(
                        make_pheno_cov,
                        gene,
                        str(tmp_path),
                        str(sample_covs_file),
                        str(celltype_covs_file),
                        str(pheno_cov_filename),
                    )
                    manage_concurrency(pheno_cov_job)
                    logging.info(f'pheno cov job for {gene}, {celltype} scheduled')

                # make cis window file
                gene_cis_filename = to_path(
                    output_path(
                        f'cis_window_files/{chromosome}/{gene_name}_{cis_window_size}bp.tsv'
                    )
                )
                if not gene_cis_filename.exists():
                    gene_cis_job = get_batch().new_python_job(
                        name=f'gene cis file: {gene}'
                    )
                    gene_cis_job.cpu(cis_job_cpu)

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

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
