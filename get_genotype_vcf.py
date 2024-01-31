#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script will

- extract common variants
- export as VCF (one per chrom)

this will be used as input for the
SAIGE-QTL association pipeline.
"""

# python modules
from cpg_utils.hail_batch import (
    dataset_path,
    get_batch,
    get_config,
    init_batch,
    output_path,
)

import click
import hail as hl
from hail.methods import export_vcf


BIOHEART_TOB_VDS = dataset_path('vds/1-0.vds')
# HAIL_IMAGE = get_config()['images']['scanpy']
BCFTOOLS_IMAGE = get_config()['images']['bcftools']

# inputs:
@click.option('--chromosomes', help=' eg chr22')
@click.command()
def main(chromosomes):
    """
    Run associaTR processing pipeline
    """

    init_batch()

    vds_path = BIOHEART_TOB_VDS
    vds = hl.vds.read_vds(vds_path)

    for chromosome in chromosomes:

        # consider only relevant chromosome
        vds = hl.vds.filter_chromosomes(vds, keep=chromosome)

        # split multiallelic loci
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)

        # densify to matrix table object
        mt = hl.vds.to_dense_mt(vds)

        # filter out loci & variant QC
        mt = mt.filter_rows(
            ~(mt.was_split)  # biallelic (exclude multiallelic)
            & (hl.len(mt.alleles) == 2)  # remove hom-ref
            & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
        )
        mt = hl.variant_qc(mt)

        # common variants only
        cv_mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.05)

        # remove fields not in the VCF
        cv_mt = cv_mt.drop('gvcf_info')

        # export to plink common variants only
        cv_vcf_path = output_path(f'{chromosome}_common_variants.vcf.bgz')
        export_vcf(cv_mt, cv_vcf_path)

        # figure out how to create index file (.csi)
        vcf_input = get_batch().read_input(cv_vcf_path)
        bcftools_job = get_batch().new_job(name=f'index vcf')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.cpu(4)
        bcftools_job.storage('15G')
        bcftools_job.command(f"bcftools index -c {vcf_input} -o {cv_vcf_path}.csi")

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
