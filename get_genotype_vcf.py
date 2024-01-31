#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script will

- extract common variants
- export as VCF (one per chrom)

this will be used as input for the
SAIGE-QTL association pipeline.


To run:

analysis-runner \
    --description "get common variant VCF" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/genotypes/" \
    python3 get_genotype_vcf.py --vds-name vds/1-0.vds chr1 chr2 chr22

"""

# python modules
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, init_batch, output_path

import click
import hail as hl
from hail.methods import export_vcf


BCFTOOLS_IMAGE = get_config()['images']['bcftools']

# inputs:
@click.option('--vds-name', help='e.g., vds/1-0.vds')
@click.argument('--chromosomes', help='e.g., chr22 chrX ', nargs=-1)
@click.command()
def main(vds_name, chromosomes):
    """
    Write genotypes as VCF
    """

    init_batch()

    vds_path = dataset_path(vds_name)
    vds = hl.vds.read_vds(vds_path)

    for chromosome in chromosomes:

        # consider only relevant chromosome
        chrom_vds = hl.vds.filter_chromosomes(vds, keep=chromosome)

        # split multiallelic loci
        chrom_vds = hl.vds.split_multi(chrom_vds, filter_changed_loci=True)

        # densify to matrix table object
        mt = hl.vds.to_dense_mt(chrom_vds)

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

        # add index file (.csi) using bcftools
        vcf_input = get_batch().read_input(cv_vcf_path)
        bcftools_job = get_batch().new_job(name='index vcf')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.cpu(4)
        bcftools_job.storage('15G')
        bcftools_job.command(f"bcftools index -c {vcf_input} -o {bcftools_job.csi}")
        get_batch().write_output(bcftools_job.csi, f'{cv_vcf_path}.csi')
    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
