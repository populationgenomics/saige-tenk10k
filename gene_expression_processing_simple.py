#!/usr/bin/env python3


"""
Hail Batch workflow to create gene expression files.
This script will:

- select genes to test based on expression (per cell type)
- build chromosome & cell type specific phenotype covariate files
- use gene info to create cis-window files
- export pheno_cov files as tsv files
- export cis_window files as tsv files

More details in README
output files in tob_wgs_genetics/saige_qtl/input

analysis-runner \
    --dataset tob-wgs \
    --access-level test \
    --output-dir 'tob_wgs_genetics/saige_qtl/hope-test-input' \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:0.0.3 \
    --description 'scRNA-seq processing batch job test' \
    python3 gene_expression_processing_simple.py \
    --celltypes=B_IN --chromosomes=chr22 \
    --gene-info-tsv=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv \
    --sample-mapping-file-path=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv \
    --expression-files-prefix=hope-test

"""

import math

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch

config = get_config()


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--gene-info-tsv')
@click.option('--expression-files-prefix')
@click.option('--sample-mapping-file-path')
@click.option('--min-pct-expr', type=int, default =5)
@click.option('--cis-window-size', type=int,default=100000)

def main(
    celltypes: str,
    chromosomes: str,
    gene_info_tsv: str,
    expression_files_prefix: str,
    sample_mapping_file_path: str,
):
    """
    Run expression processing pipeline
    """

    b = get_batch()
    # create phenotype covariate files
    smf_df = pd.read_csv(sample_mapping_file_path, sep='\t')
    gene_info_df = pd.read_csv(gene_info_tsv, sep='\t')

    
    # set jobs running
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
