#!/usr/bin/env python3


"""
Hail Batch workflow to get_chrom_celltype_expression_and_filter.
This script will:

- select genes to test based on expression (per cell type) and filter out lowly expressed genes 

More details in README
output files in tob_wgs_genetics/saige_qtl/input

   

    analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files" \
    --access-level "test" \
    --output-dir "tob_wgs_genetics/saige_qtl/hope-test-input" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
     gene_expression_processing_part1.py  --celltypes=B_IN --chromosomes=chr22 \
    --gene-info-tsv=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv \
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


from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

from cpg_utils.hail_batch import copy_common_env, dataset_path, output_path

config = get_config()


SCANPY_IMAGE = config['images']['scanpy']


def get_chrom_celltype_expression_and_filter(
    gene_info_df,
    expression_files_prefix: str,  # tob_wgs_genetics/saige_qtl/input/
    chromosome: str,
    cell_type: str,
    min_pct: int,
    h5ad_ofile_path: str
):
    """Extracts relevant expression info AND remove genes with low expression across cells 

    Input:
    - chromosome & cell type of interest
    - path to (single-cell) expression files, one tsv file per cell type,
    rows are cells and columns are genes
    - path to dataframe containing gene info, for each gene (=row),
    specifies chrom, start, end and strand

    Output: GCS path to FILTERED expression adata object for only relevant genes
    """

    # first line is where the file is now,
    # but (second line) the files will eventually be in the below folder
    # and split by cell type (at least this all naive B cells only)
    expression_h5ad_path = to_path(
        dataset_path(
            f'scrna-seq/CellRegMap_input_files/expression_objects/sce22.h5ad'
        )
    ).copy('here.h5ad')
    expression_adata = scanpy.read(expression_h5ad_path)

    # select only genes on relevant chromosome
    genes_chrom = gene_info_df[gene_info_df['chr'] == chromosome].gene_name
    # return expression for the correct chromosomes only
    expression_adata = expression_adata[:, expression_adata.var_names.isin(genes_chrom)]
    
    #filter lowly expressed genes 
    n_all_cells = len(expression_adata.obs.index)
    min_cells_input = math.ceil((n_all_cells * min_pct) / 100)
    scanpy.pp.filter_genes(expression_adata, min_cells=min_cells_input)
    assert isinstance(expression_adata, scanpy.AnnData)
    print(expression_adata)

    # write expression_adata to tmp file path
    expression_adata.write_h5ad(str(h5ad_ofile_path))

    # write genes array to GCS directly (because we can!)
    with to_path(output_path(f'{chromosome}_{cell_type}_filtered_genes.json')).open('w') as write_handle:
        json.dump(list(expression_adata.var_names), write_handle)



@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--gene-info-tsv')
@click.option('--expression-files-prefix')
@click.option('--min-pct-expr', type=int, default =5)

def main(
    celltypes: str,
    chromosomes: str,
    gene_info_tsv: str,
    expression_files_prefix: str,
    min_pct_expr: int,
):
    """
    Run expression processing pipeline
    """
    config = get_config()
    b = get_batch()

    logging.info(f'Cell types to run: {celltypes}')
    logging.info(f'Chromosomes to run: {chromosomes}')

    gene_info_df = pd.read_csv(gene_info_tsv, sep='\t')

    for celltype in celltypes.split(','):
        
        for chromosome in chromosomes.split(','):
            j = b.new_python_job(name=f'Get expression (celltype:{celltype} and chromosome:{chromosome}), then filter lowly exp. genes')
            j.storage('20G')
            j.cpu(8)
            j.image(config['workflow']['driver_image'])
            j.call(get_chrom_celltype_expression_and_filter,gene_info_df,expression_files_prefix,chromosome,celltype,min_pct_expr,j.ofile)
            j.ofile.add_extension('.h5ad')
            b.write_output(j.ofile, output_path(f'filtered_{celltype}_{chromosome}.h5ad'))

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter