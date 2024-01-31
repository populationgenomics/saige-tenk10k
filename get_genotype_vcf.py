#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script will

- extract common variants
- export as VCF (one per chrom)
- also create a plink file for a subset of variants
for variance ratio estimation

this will be used as input for the
SAIGE-QTL association pipeline.


To run:

analysis-runner \
    --description "get common variant VCF" \
    --dataset "bioheart" \
    --access-level "test" \
    --output-dir "saige-qtl/input_files/genotypes/" \
    python3 get_genotype_vcf.py --vds-version 1-0 --chromosomes chr1,chr2,chr22

"""

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, init_batch, output_path

import click
import random
import hail as hl
from hail.methods import export_plink, export_vcf


BCFTOOLS_IMAGE = get_config()['images']['bcftools']


def can_reuse(path: str):
    if path and to_path(path).exists():
        return True
    return False


# inputs:
@click.option('--vds-version', help=' e.g., 1-0 ')
@click.option('--chromosomes', help=' e.g., chr22,chrX ')
@click.option('--vre-mac-threshold', default=20)
@click.option('--vre-n-markers', default=2000)
@click.command()
def main(vds_version, chromosomes, vre_mac_threshold, vre_n_markers):
    """
    Write genotypes as VCF
    """

    init_batch()

    vds_path = dataset_path(f'vds/{vds_version}.vds')
    vds = hl.vds.read_vds(vds_path)

    for chromosome in chromosomes.split(','):

        # create path and check if it exists already
        cv_vcf_path = output_path(
            f'vds{vds_version}/{chromosome}_common_variants.vcf.bgz'
        )
        if not can_reuse(cv_vcf_path):

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
            export_vcf(cv_mt, cv_vcf_path)

            # add index file (.csi) using bcftools
            vcf_input = get_batch().read_input(cv_vcf_path)
            bcftools_job = get_batch().new_job(name='index vcf')
            bcftools_job.image(BCFTOOLS_IMAGE)
            bcftools_job.cpu(4)
            bcftools_job.storage('15G')
            bcftools_job.command(f"bcftools index -c {vcf_input} -o {bcftools_job.csi}")
            get_batch().write_output(bcftools_job.csi, f'{cv_vcf_path}.csi')

    # subset variants for variance ratio estimation
    vre_plink_path = output_path(f'vds{vds_version}/vre_plink_2000_variants')
    if not can_reuse(vre_plink_path):
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)
        mt = hl.vds.to_dense_mt(vds)
        mt = mt.filter_rows(
            ~(mt.was_split)  # biallelic (exclude multiallelic)
            & (hl.len(mt.alleles) == 2)  # remove hom-ref
            & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs (exclude indels)
        )
        mt = hl.variant_qc(mt)
        # minor allele count (MAC) > {vre_mac_threshold}
        vre_mt = mt.filter_rows(mt.variant_qc.AC[0] > {vre_mac_threshold})
        n_ac_vars = vre_mt.count()[0]  # to avoid evaluating this 2X
        if n_ac_vars == 0:
            return
        # since pruning is very costly, subset first a bit
        random.seed(0)
        vre_mt = vre_mt.sample_rows(p=0.01)
        # perform LD pruning
        pruned_variant_table = hl.ld_prune(vre_mt.GT, r2=0.2, bp_window_size=500000)
        vre_mt = vre_mt.filter_rows(hl.is_defined(pruned_variant_table[vre_mt.row_key]))
        # randomly sample {vre_n_markers} variants
        random.seed(0)
        vre_mt = vre_mt.sample_rows(({vre_n_markers} * 1.1) / vre_mt.count()[0])
        vre_mt = vre_mt.head({vre_n_markers})

        # export to plink common variants only for sparse GRM
        export_plink(vre_mt, vre_plink_path, ind_id=vre_mt.s)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
