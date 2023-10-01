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
    gene_expression_processing_simple.py \

"""

from cpg_workflows.batch import get_batch



@click.command()

def main():
    b = get_batch()

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
