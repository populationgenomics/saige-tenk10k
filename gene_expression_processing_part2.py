#!/usr/bin/env python3


"""
Hail Batch workflow to create gene expression files.
This script will:

- build chromosome & cell type specific phenotype covariate files
- use gene info to create cis-window files
- export pheno_cov files as tsv files
- export cis_window files as tsv files

More details in README
output files in tob_wgs_genetics/saige_qtl/input

   

    analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files" \
    --access-level "test" \
    --output-dir "tob_wgs_genetics/saige_qtl/hope-test-input" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
     gene_expression_processing_part2.py  --celltypes=B_IN --chromosomes=chr22 \
    --gene-info-tsv=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv \
    --sample-mapping-file-path=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv \
    --expression-files-prefix=hope-test

"""

import os
import logging
import math
import sys
#from tkinter.tix import CELL

import click
import hail as hl
import pandas as pd
from cpg_workflows.batch import get_batch
import hailtop.batch as hb
import scanpy
import json
import numpy as np


from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

from cpg_utils.hail_batch import copy_common_env, dataset_path, output_path

config = get_config()


SCANPY_IMAGE = config['images']['scanpy']

def get_celltype_covariates(
    expression_files_prefix: str,
    cell_type: str,
):
    """Obtain cell type specific covariates

    Input:
    - cell type of interest
    - covariate files prefix

    Output: covariate df for cell type of interest
    """
    covs_tsv_path = dataset_path('tob_wgs_genetics/saige_qtl/input/covariate_chr22_B_IN_tester_final.csv')
    covs_df = pd.read_csv(covs_tsv_path, sep=',', index_col=0)
    return covs_df


def build_pheno_cov_filename(
    gene_name, expression_files_prefix, celltype,chromosome, sample_mapping_file_path, ofile_path: str  # sample mapping file
):
    """
    Combine files to build final input
    reduce to one file per gene

    Input:
    - Expression anndata (filtered)
    - Sample mapping file path

    """
    #read in filtered anndata file: (test by taking it out of the for loop)
    filtered_h5ad_path = to_path((output_path(f'filtered_{celltype}_{chromosome}.h5ad'))).copy('here.h5ad')
    expression_adata = scanpy.read(filtered_h5ad_path)

    smf_df = pd.read_csv(sample_mapping_file_path, sep='\t')
    cov_df = get_celltype_covariates(expression_files_prefix, celltype)

    gene_index = np.where(expression_adata.raw.var.index == gene_name)[0][0]
    gene_mat = expression_adata.raw.X[:,gene_index].todense()

    expression_df = pd.DataFrame(
    data=gene_mat,
    index=expression_adata.obs.index,  #cell IDs
    columns=expression_adata.raw.var.index[np.where(expression_adata.raw.var.index == gene_name)]
)
    expression_df['OneK1K_ID'] = expression_df.obs['individual']

    # Reset the index and make it a column
    expression_df.reset_index(inplace=True)

    # Rename the columns 
    expression_df.columns = ['Cell_ID', f'{gene_name}_raw_count', 'OneK1K_ID']

    expression_df = pd.merge(expression_df, smf_df, on='OneK1K_ID', how='left')

    pheno_cov_df = pd.merge(expression_df, cov_df, on = "Cell_ID", how = "inner")

    pheno_cov_df.to_csv(str(ofile_path), index=False)

def get_gene_cis_file(gene_info_df, gene: str, window_size: int):
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
@click.option('--gene-info-tsv')
@click.option('--expression-files-prefix')
@click.option('--sample-mapping-file-path')
@click.option('--min-pct-expr', type=int, default =5)
@click.option('--cis-window-size', type=int,default=100000)
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
    gene_info_tsv: str,
    expression_files_prefix: str,
    sample_mapping_file_path: str,
    min_pct_expr: int,
    cis_window_size: int,
    max_gene_concurrency=int
):
    """
    Run expression processing pipeline
    """
    config = get_config()
    b = get_batch()
    

    logging.info(f'Cell types to run: {celltypes}')
    logging.info(f'Chromosomes to run: {chromosomes}')

    # create phenotype covariate files
    #gene_info_df = pd.read_csv(gene_info_tsv, sep='\t')

    for celltype in celltypes.split(','):
        for chromosome in chromosomes.split(','):
            #load array containing the genes from filtered anndata 
            with to_path(output_path(f'{chromosome}_{celltype}_filtered_genes.json')).open('r') as read_handle:
                genes = json.load(read_handle)

            # combine files for each gene
            # pylint: disable=no-member
            for gene in genes:
                pheno_cov_job = b.new_python_job(name=f'Build phenotype-covariate files for {gene} [{celltype};{chromosome}]')
                pheno_cov_job.storage('35G')
                pheno_cov_job.cpu(8)
                pheno_cov_job.image(config['workflow']['driver_image'])
                pheno_cov_job.call(
                    build_pheno_cov_filename,gene,"", celltype, chromosome, sample_mapping_file_path,pheno_cov_job.ofile)
                pheno_cov_job.ofile.add_extension('.csv')
                b.write_output(pheno_cov_job.ofile, output_path(
                      f'input_files/pheno_cov_files/{gene}_{celltype}.csv'
                    ))

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
