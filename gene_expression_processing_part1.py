#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter,import-error, line-too-long

"""
PART 1 of 2 of the 'Gene expression processing workflow' 

Aims: 
- Get chromosome-level and cell-type specific expression data
- and filter lowly expressed genes.

More details in README
output files in tob_wgs_genetics/saige_qtl/input

    analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files" \
    --access-level "test" \
    --output-dir "tob_wgs_genetics/saige_qtl/hope-test-input" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
     gene_expression_processing_part1.py  \
    --celltypes=B_IN,CD4_NC,CD4_ET,CD4_SOX4,CD8_ET,CD8_NC,CD8_S100B,\
    NK,NK_R,Plasma,B_MEM,B_IN,MonoC,MonoNC,DC \
    --chromosomes=chr22 \
    --gene-annotation-file=gs://cpg-tob-wgs-main/tob_wgs_genetics/gencode.v42.annotation.gff3.gz \
    --expression-files-prefix=scrna-seq/CellRegMap_input_files/expression_objects

"""

import logging
import math
import json

import click
import pandas as pd
from cpg_workflows.batch import get_batch
import scanpy


from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, output_path

config = get_config()

SCANPY_IMAGE = config['images']['scanpy']


def gene_info(x):
    """Helper function to extract ENSG and gene_name from a GFF3 annotation file"""
    g_name = list(filter(lambda x: 'gene_name' in x, x.split(';')))[0].split('=')[1]
    g_id = list(filter(lambda x: 'gene_id' in x, x.split(';')))[0].split('=')[1]
    g_id = g_id.split('.')[0]  # removes the version number from ENSG ids
    return (g_name, g_id)


def get_chrom_celltype_expression_and_filter(
    gene_annotation_file: str,
    expression_files_prefix: str,
    chromosome: str,
    cell_type: str,
    min_pct: int,
    h5ad_ofile_path: str,
):
    """Extracts relevant expression info AND remove genes with low expression across cells

    Input:
    - path to GENCODE annotation file in GFF3 format
    - path to (single-cell) expression files (anndata objects, separated by chromosome)
    - chromosome & cell type of interest
    - path to (single-cell) expression files, one tsv file per cell type,
    rows are cells and columns are genes


    Output: GCS path to FILTERED expression adata object for only relevant genes
    """
    chromosome_number = chromosome[3:]

    expression_h5ad_path = to_path(
        dataset_path(f'{expression_files_prefix}/sce{chromosome_number}.h5ad')
    ).copy('here.h5ad')
    expression_adata = scanpy.read(expression_h5ad_path)

    # Manual editing of a mislabeled individual ID
    condition = (expression_adata.obs['individual'] == '870_871') & (
        expression_adata.obs['latent'] == 0
    )
    # Update the 'individual' attribute for cells that meet the condition
    expression_adata.obs.loc[condition, 'individual'] = '966_967'

    # cell_type - labels to number mapping
    cell_type_mapping = {
        'CD4_NC': 0,
        'CD4_ET': 1,
        'CD4_SOX4': 2,
        'CD8_ET': 3,
        'CD8_NC': 4,
        'CD8_S100B': 5,
        'DC': 6,
        'Plasma': 8,
        'MonoC': 9,
        'MonoNC': 10,
        'B_MEM': 12,
        'B_IN': 13,
        'NK': 14,
        'NK_R': 15,
    }

    # filter expression_adata based on cell type
    cell_type_numerical_label = cell_type_mapping[cell_type]
    expression_adata = expression_adata[
        expression_adata.obs.cell_type == cell_type_numerical_label
    ]

    # Reads and wrangles GENCODE annotation to extract ENSG ID and gene name mappings
    gencode = pd.read_table(
        gene_annotation_file,
        comment='#',
        sep='\t',
        names=[
            'seqname',
            'source',
            'feature',
            'start',
            'end',
            'score',
            'strand',
            'frame',
            'attribute',
        ],
    )
    gencode_genes = (
        gencode[(gencode.feature == 'gene')][['seqname', 'start', 'end', 'attribute']]
        .copy()
        .reset_index()
        .drop('index', axis=1)
    )  # subsets for gene annotations
    gencode_genes['gene_name'], gencode_genes['ENSG'] = zip(
        *gencode_genes.attribute.apply(gene_info)
    )

    # subset gencode annotation file for relevant chromosome
    gencode_genes = gencode_genes[gencode_genes['seqname'] == chromosome]

    # Intersects gene annotations in single cell dataset with gencode gene annotations
    expression_adata = expression_adata[
        :, expression_adata.var_names.isin(gencode_genes['gene_name'])
    ]

    # filter lowly expressed genes
    n_all_cells = len(expression_adata.obs.index)
    min_cells_input = math.ceil((n_all_cells * min_pct) / 100)
    scanpy.pp.filter_genes(expression_adata, min_cells=min_cells_input)
    assert isinstance(expression_adata, scanpy.AnnData)
    print(expression_adata)

    # write expression_adata to tmp file path
    expression_adata.write_h5ad(str(h5ad_ofile_path))

    # write genes array to GCS directly
    with to_path(
        output_path(f'{cell_type}/{chromosome}_{cell_type}_filtered_genes.json')
    ).open('w') as write_handle:
        json.dump(list(expression_adata.var_names), write_handle)


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes', help='example chr22')
@click.option('--gene-annotation-file', help='Gencode annotation file in GFF3 format')
@click.option('--expression-files-prefix')
@click.option('--min-pct-expr', type=int, default=5)
def main(
    celltypes: str,
    chromosomes: str,
    gene_annotation_file: str,
    expression_files_prefix: str,
    min_pct_expr: int,
):
    """Extracts relevant expression info AND remove genes with low expression across by cell type and chromosome"""

    b = get_batch()

    logging.info(f'Cell types to run: {celltypes}')
    logging.info(f'Chromosomes to run: {chromosomes}')

    for celltype in celltypes.split(','):

        for chromosome in chromosomes.split(','):
            j = b.new_python_job(
                name=f'Get expression (celltype:{celltype}; chromosome:{chromosome}), then filter lowly exp. genes'
            )
            j.storage('4G')
            j.cpu(4)
            j.image(config['workflow']['driver_image'])
            j.call(
                get_chrom_celltype_expression_and_filter,
                gene_annotation_file,
                expression_files_prefix,
                chromosome,
                celltype,
                min_pct_expr,
                j.ofile,
            )
            j.ofile.add_extension('.h5ad')
            b.write_output(
                j.ofile, output_path(f'filtered_{celltype}_{chromosome}.h5ad')
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()
