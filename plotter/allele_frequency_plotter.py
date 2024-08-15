#!/usr/bin/env python3
"""
This script plots a histogram of allele frequencies for a given dataset

To run:

analysis-runner \
    --description "plot frequencies histogram" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/" \
    python3 allele_frequency_plotter.py \
        --vds-path=gs://cpg-bioheart-test/vds/bioheart1-0.vds \
        --title='frequencies BioHEART'

vds path tob: gs://cpg-tob-wgs-test/vds/tob-wgs1-0.vds
vds path bioheart: gs://cpg-bioheart-test/vds/bioheart1-0.vds
vds path tob+bioheart: gs://cpg-bioheart-test/vds/tenk10k1-0.vds
"""

import click

import hail as hl

from cpg_utils.hail_batch import init_batch, output_path
from bokeh.plotting import output_file, save


@click.command()
@click.option('--vds-path', required=True)
@click.option('--title', default='Minor allele frequency')
@click.option('--exclude-multiallelic', default=False)
@click.option('--exclude-indels', default=False)
def plot_frequencies(
    vds_path: str,
    title: str,
    exclude_multiallelic: bool,
    exclude_indels: bool,
):
    """
    reads the VDS, converts to MT,
    if set to true exclude indels and multiallelic snps
    and plots the frequency distributions
    of the remaining variants
    """
    # read VDS object (WGS data)
    init_batch()
    vds_name = vds_path.split('/')[-1].split('.')[0]
    print(f'Plotting variants frequencies for {vds_name}')
    vds = hl.vds.read_vds(vds_path)
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    # convert to hail matrix table
    mt = hl.vds.to_dense_mt(vds)

    # filter out loci & variant QC
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # remove hom-ref
    if exclude_multiallelic:  # biallelic only (exclude multiallelic)
        print('Excluding multiallelic loci')
        mt = mt.filter_rows(~(mt.was_split))
    if exclude_indels:  # SNPs only (exclude indels)
        print('Excluding indels')
        mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

    # compute allele frequencies as part of variant qc
    mt = hl.variant_qc(mt)

    p = hl.plot.histogram(mt.variant_qc.AF[1], legend=title)
    output_file('local_histo.html')
    save(p)
    gcs_path_p = output_path(f'plots/frequency_histo/{vds_name}.html', 'analysis')
    hl.hadoop_copy('local_histo.html', gcs_path_p)


if __name__ == '__main__':
    plot_frequencies()
