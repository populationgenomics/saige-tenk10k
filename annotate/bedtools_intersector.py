#!/usr/bin/env python3


"""
Intersect with Arthur's annotations
analysis-runner --dataset "bioheart" \
    --description "bedtools_intersector" \
    --access-level "test" \
    --output-dir "saige-qtl/arthur" \
        annotate/bedtools_intersector.py --catalog-path='gs://cpg-bioheart-test/str/ncAnnot.v0.14.jul2024.bed' \
        --query-path='gs://cpg-bioheart-test/saige-qtl/arthur/SH2D2A_cis_region.tsv'


essentially a copy of
https://github.com/populationgenomics/sv-workflows/blob/57af2557cca1be418428351b42210639ecaa70f0/str/annotate/bedtools_intersector.py

"""
import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch


@click.command()
@click.option('--query-path', required=True)
@click.option(
    '--catalog-path', default='gs://cpg-bioheart-test/str/ncAnnot.v0.14.jul2024.bed'
)
def bedtools_intersect(
    query_path: str,
    catalog_path: str,
):
    # Initializing Batch
    b = get_batch('Bedtools intersect')

    # provision a new job
    bedtools_job = b.new_job(name='Bedtools intersect')
    bedtools_job.image(get_config()['images']['bedtools'])
    bedtools_job.storage('50G')
    bedtools_job.memory('highmem')
    bedtools_job.cpu(16)

    # read input files
    catalog = b.read_input(catalog_path)
    tr = b.read_input(query_path)

    # set the job command
    bedtools_job.command(
        f'bedtools intersect -a {tr} -b {catalog} -wo > {bedtools_job.ofile}'
    )

    # write output to GCP bucket for this dataset
    output_file = query_path.replace('.tsv', '_arthur_intersect.bed')
    b.write_output(bedtools_job.ofile, output_file)

    b.run(wait=False)


if __name__ == '__main__':
    bedtools_intersect()
