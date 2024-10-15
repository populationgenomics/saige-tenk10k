#!/usr/bin/env python3

"""
Reads in a single VDS and generates per-chromosome MTs

analysis-runner \
   --description "generate chr subsets" \
   --dataset "bioheart" \
   --access-level "standard" \
   --output-dir "generate_chr_subsets" \
   python3 generate_chrom_vds_subsets.py -i <VDS> -o gs://cpg-bucket/path/to/subsets --chroms chr2 chr3 chr4
"""

import logging
from argparse import ArgumentParser
from os.path import join

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch


def make_chrom_subset(vds_path: str, chrom: str, out_path: str):
    """
    Make single CHR subset
    """
    import random
    import time

    from cpg_utils.hail_batch import init_batch

    time.sleep(random.randint(0, 400))

    init_batch()

    # restrict to one chromosome, and convert to MT
    vds = hl.vds.read_vds(vds_path)
    chrom_vds = hl.vds.filter_chromosomes(vds, keep=chrom)
    chrom_mt = hl.vds.to_dense_mt(chrom_vds)
    chrom_mt.write(out_path)


def main():
    """
    shard the input VDS into per-chromosome MTs
    """

    parser = ArgumentParser()
    parser.add_argument('-i', required=True, help='input VDS')
    parser.add_argument('-o', required=True, help='output MT root (folder)')
    parser.add_argument('--chroms', nargs='+', default=None)
    args = parser.parse_args()

    init_batch()

    chroms = args.chroms if args.chroms else [f'chr{chr}' for chr in range(1, 23)]

    # loop over chromosomes
    for chrom in chroms:

        out_path = join(args.o, f'{chrom}_subset.mt')

        # check for successful prior write
        success_path = join(out_path, '_SUCCESS')

        if to_path(success_path).exists():
            logging.info(f'Skipping {chrom}, as {success_path} exists')

        chrom_job = get_batch().new_python_job(name=f'Make subset MT for {chrom}')
        chrom_job.storage('50Gi')
        chrom_job.call(
            make_chrom_subset, vds_path=args.i, chrom=chrom, out_path=out_path
        )

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
