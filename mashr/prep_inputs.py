#!/usr/bin/env python3

"""
This script prepares eSTR inputs from every cell type.

analysis-runner --dataset "bioheart" \
    --description "Prepare inputs for mashr" \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --access-level "test" \
    --output-dir "saige-qtl/mashr" \
    prep_inputs.py

copy of https://github.com/populationgenomics/sv-workflows/blob/mashR/str/mashR/prep_inputs.py
"""

import os
import sys

import numpy as np
import pandas as pd

from cpg_utils.hail_batch import get_batch
from cpg_utils import to_path


def celltype_chrom_parser(celltype, chrom, estrs_coord_chrom):
    celltype_df = pd.DataFrame()
    for index, row in estrs_coord_chrom.iterrows():
        chrom = row['chr']
        gene_name = row['gene_name']
        pos = row['pos']
        motif = row['motif']
        ref_len = row['ref_len']
        try:
            df = pd.read_csv(
                f'gs://cpg-bioheart-main-analysis/str/associatr/common_variants_snps/tob_n1055_and_bioheart_n990/meta_results/meta_results/{celltype}/{chrom}/{gene_name}_100000bp_meta_results.tsv',
                sep='\t',
            )
            df = df[
                (df['pos'] == pos) & (df['motif'] == motif) & (df['ref_len'] == ref_len)
            ]
            beta = df['coeff_meta'].iloc[0]
            se = df['se_meta'].iloc[0]
            # save the beta and se into a row:
            beta_se = pd.DataFrame(
                {
                    'chrom': chrom,
                    'pos': pos,
                    'motif': motif,
                    'ref_len': ref_len,
                    'gene': gene_name,
                    f'{celltype}_beta': beta,
                    f'{celltype}_se': se,
                },
                index=[0],
            )
            celltype_df = pd.concat([celltype_df, beta_se], axis=0)
        except:
            continue

    celltype_df.to_csv(
        f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/estrs_beta_se/{celltype}/{chrom}/beta_se.tsv',
        sep='\t',
        index=False,
    )


def celltype_chrom_parser_null(celltype, chrom):
    gene_files = list(
        to_path(
            f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/{celltype}/chr{chrom}'
        ).rglob('*.tsv')
    )
    master_df = pd.DataFrame()
    for gene_file in gene_files:
        gene_name = str(gene_file).split('/')[-1].split('_')[0]
        df = pd.read_csv(gene_file, sep='\t')
        df['gene'] = gene_name
        df[f'{celltype}_beta'] = df['coeff_meta']
        df[f'{celltype}_se'] = df['se_meta']
        df = df[
            [
                'chr',
                'pos',
                'motif',
                'ref_len',
                'gene',
                f'{celltype}_beta',
                f'{celltype}_se',
            ]
        ]
        master_df = pd.concat([master_df, df], axis=0)

    master_df.to_csv(
        f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/chr22_null_beta_se/{celltype}/chr{chrom}/beta_se.tsv',
        sep='\t',
        index=False,
    )


def main():
    # b = get_batch(name='Prep eSTRs for mashr')
    cell_types = 'CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC'

    celltypes = cell_types.split(',')
    # load in the list of eSTRs passing FDR <5% across all cell types:
    # estrs_coord = pd.read_csv(
    #    'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/estrs_coord_gene.csv',
    # )
    master_df = pd.DataFrame()
    for celltype in celltypes:

        # for chrom in [22]:
        # df = pd.read_csv(
        #    f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/estrs_beta_se/{cell}/chr{chrom}/beta_se.tsv',
        #    sep='\t',
        # )
        df = pd.read_csv(
            f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/chr22_null_beta_se/{celltype}/chr22/beta_se.tsv',
            sep='\t',
        )
        if master_df.empty:
            master_df = df
        else:
            master_df = master_df.merge(
                df, on=['chr', 'pos', 'motif', 'ref_len', 'gene'], how='inner'
            )
            # estrs_coord_chrom = estrs_coord[estrs_coord['chr'] == f'chr{chrom}']
            # if to_path(f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/estrs_beta_se/{cell}/{chrom}/beta_se.tsv').exists():
            #    continue
            # job = b.new_python_job(f'Prep eSTRs for mashr {cell} {chrom}')
            # job.cpu(0.25)
            # job.call(cell_chrom_parser_null, cell, chrom)
    master_df.to_csv(
        f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/chr22_nullbeta_se/chr22/all_cell_chr22_beta_se.tsv',
        sep='\t',
        index=False,
    )
    # b.run(wait=False)


if __name__ == '__main__':
    main()
