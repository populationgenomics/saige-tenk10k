#!/usr/bin/env python3

"""
Copy of https://github.com/populationgenomics/sv-workflows/blob/main/str/coloc/coloc_runner.py
to rework for common variant SAIGE-QTL results

This script performs SNP-only colocalisation analysis betweeen eGenes identified by single-cell eQTL analysis (using SAIGE-QTL) and GWAS signals.
Assumes that the SNP GWAS data has been pre-processed with the following columns: 'chromosome', 'position' (hg38 bp), 'snp'(chromosome_position_refallele_effectallele), 'beta', 'varbeta'

1) Identify eGenes where at least one eQTL has pval < 5e-8
2) Extract the SNP GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene (if the SNP GWAS data has at least one variant with pval <5e-8)
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by SAIGE-QTL analysis" \
    --access-level "test" \
    --memory='8G' \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/" \
    coloc/coloc_runner.py \
    --snp-gwas-file=gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/ibd_EAS_EUR_SiKJEF_meta_IBD.tsv \
    --pheno-output-name="ibd_liu2023" \
    --celltypes "NK"

"""

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, image_path, output_path


def coloc_runner(gwas, eqtl_file_path, celltype, pheno_output_name):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(coloc)')
    ro.r('library(tidyverse)')

    with (ro.default_converter + pandas2ri.converter).context():
        gwas_r = ro.conversion.get_conversion().py2rpy(gwas)
    ro.globalenv['gwas_r'] = gwas_r
    ro.r(
        '''
    gwas_r = gwas_r %>% select(beta, varbeta, position,snp)
    gwas_r = gwas_r %>% filter((beta!=0) | (varbeta!=0))
    gwas_r = gwas_r %>% distinct(snp, .keep_all = TRUE)
    gwas_r = gwas_r%>% as.list()
    gwas_r$type = 'cc'

    ''',
    )
    # TODO: update eQTL df parsing to SAIGE-QTL results
    eqtl = pd.read_csv(
        eqtl_file_path,
        sep='\t',
    )
    eqtl['beta'] = eqtl['BETA']
    eqtl['se'] = eqtl['SE']
    eqtl['position'] = eqtl['POS']
    eqtl['snp'] = eqtl['MarkerID'].str.replace(':', '_', regex=False)
    # eqtl['beta'] = eqtl['coeff_meta']
    # eqtl['se'] = eqtl['se_meta']
    # eqtl['position'] = eqtl['pos']
    # eqtl['snp'] = eqtl['chr'] + '_' + eqtl['position'].astype(str) + '_' + eqtl['motif']
    # eqtl['snp'] = eqtl['snp'].str.replace('-', '_', regex=False)
    gene = eqtl_file_path.split('/')[-1].replace(f'{celltype}_', '').replace('_cis', '')
    # gene = eqtl_file_path.split('/')[-1].split('_')[0]
    with (ro.default_converter + pandas2ri.converter).context():
        eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
    ro.globalenv['eqtl_r'] = eqtl_r
    ro.globalenv['gene'] = gene
    ro.r(
        '''
    eqtl_r = eqtl_r %>% filter(!is.na(beta))
    eqtl_r = eqtl_r %>% distinct(snp, .keep_all = TRUE)
    eqtl_r$varbeta = eqtl_r$se**2
    eqtl_r$position = eqtl_r$pos
    eqtl_r = eqtl_r %>% select(beta, varbeta, position, snp)

    eqtl_r = eqtl_r %>% as.list()
    eqtl_r$type = 'quant'
    eqtl_r$sdY = 1


    my.res <- coloc.abf(dataset1=gwas_r,
                    dataset2=eqtl_r)

    p_df <- data.frame(gene,my.res$summary[1], my.res$summary[2], my.res$summary[3], my.res$summary[4], my.res$summary[5], my.res$summary[6])
    names(p_df) <- c('gene', 'nsnps_coloc_tested','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_p4_df = ro.conversion.get_conversion().rpy2py(ro.r('p_df'))

    # add cell type and chrom annotation to df
    pd_p4_df['celltype'] = celltype
    pd_p4_df['chrom'] = eqtl['chr'].iloc[0]

    # write to GCS
    pd_p4_df.to_csv(
        output_path(
            f"coloc-snp-only/sig_str_filter_only/{pheno_output_name}/{celltype}/{gene}_100kb.tsv",
            'analysis',
        ),
        sep='\t',
        index=False,
    )


@click.option(
    '--egenes-file',
    help='Path to the eGenes file with FINEMAP and SUSIE probabilities',
    default='gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/all_cell_types_all_genes_sig_only.tsv',
)
@click.option(
    '--snp-cis-dir',
    help='Path to the directory containing the SNP cis results',
    default='gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files',
)
@click.option(
    '--snp-gwas-file',
    help='Path to the SNP GWAS file',
    default='gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv',
)
@click.option('--celltypes', help='Cell types to run', default='ASDC')
@click.option('--cis-window-size', help='Cis window size used', default=100000)
@click.option(
    '--max-parallel-jobs', help='Maximum number of parallel jobs to run', default=500
)
@click.option(
    '--pheno-output-name', help='Phenotype output name', default='covid_GCST011071'
)
@click.option('--job-cpu', help='Number of CPUs to use for each job', default=0.25)
@click.command()
def main(
    snp_cis_dir,
    egenes_file,
    celltypes,
    snp_gwas_file,
    pheno_output_name,
    max_parallel_jobs,
    cis_window_size,
    job_cpu,
):
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    # read in gene annotation file
    var_table = pd.read_csv(
        'gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv',
    )
    hg38_map = pd.read_csv(
        snp_gwas_file,
        sep='\t',
    )

    # read in eGenes file
    # TODO: update to take in SAIGE-QTL results
    result_df_cfm = pd.read_csv(
        egenes_file,
        sep='\t',
        usecols=[
            'chr',
            'pos',
            'pval_meta',
            'motif',
            'susie_pip',
            'gene',
            'finemap_prob',
            'celltype',
            'ref_len',
        ],
    )

    result_df_cfm['variant_type'] = (
        result_df_cfm['motif'].str.contains('-').map({True: 'SNV', False: 'STR'})
    )
    result_df_cfm_str = result_df_cfm[
        result_df_cfm['variant_type'] == 'STR'
    ]  # filter for STRs
    result_df_cfm_str = result_df_cfm_str[
        result_df_cfm_str['pval_meta'] < 5e-8
    ]  # filter for STRs with p-value < 5e-8
    result_df_cfm_str = result_df_cfm_str.drop_duplicates(
        subset=['gene', 'celltype'],
    )  # drop duplicates (ie pull out the distinct genes in each celltype)
    result_df_cfm_str['gene'] = result_df_cfm_str['gene'].str.replace(
        '.tsv',
        '',
        regex=False,
    )  # remove .tsv from gene names (artefact of the data file)
    b = get_batch(name=f'Run coloc:{pheno_output_name}')

    for celltype in celltypes.split(','):
        result_df_cfm_str_celltype = result_df_cfm_str[
            result_df_cfm_str['celltype'] == celltype
        ]  # filter for the celltype of interest
        for gene in result_df_cfm_str_celltype['gene']:
            chrom = result_df_cfm_str_celltype[
                result_df_cfm_str_celltype['gene'] == gene
            ]['chr'].iloc[0]
            if to_path(
                output_path(
                    f"coloc-snp-only/sig_str_filter_only/{pheno_output_name}/{celltype}/{gene}_{cis_window_size}.tsv",
                    'analysis',
                ),
            ).exists():
                continue

            # TO DO: figure out whether or not to include to_path here
            eqtl_results_file = (
                f'{snp_cis_dir}/{celltype}/{chrom}/{celltype}_{gene}_cis'
            )
            # if to_path(
            #     f'{snp_cis_dir}/{celltype}/{chrom}/{celltype}_{gene}_cis'
            # ).exists():
            if to_path(eqtl_results_file).exists():
                print('Cis results for ' + gene + ' exist: proceed with coloc')

                # extract the coordinates for the cis-window (gene +/- 100kB)
                gene_table = var_table[var_table['gene_ids'] == gene]
                start = int(gene_table['start'].astype(int)) - cis_window_size
                end = int(gene_table['end'].astype(int)) + cis_window_size
                chrom = gene_table['chr'].iloc[0]
                hg38_map_chr = hg38_map[hg38_map['chromosome'] == (chrom)]
                hg38_map_chr_start = hg38_map_chr[hg38_map_chr['position'] >= start]
                hg38_map_chr_start_end = hg38_map_chr_start[
                    hg38_map_chr_start['position'] <= end
                ]
                if hg38_map_chr_start_end.empty:
                    print(
                        'No SNP GWAS data for '
                        + gene
                        + ' in the cis-window: skipping....'
                    )
                    continue
                # check if the p-value column contains at least one value which is <=5e-8:
                # if hg38_map_chr_start_end['p_value'].min() > 5e-8:
                # print('No significant SNP GWAS data for ' + gene + ' in the cis-window: skipping....')
                # continue
                # print('Extracted SNP GWAS data for ' + gene)

                # run coloc
                coloc_job = b.new_python_job(
                    f'Coloc for {gene}: {celltype}',
                )
                # f'{snp_cis_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv' ???
                coloc_job.image(image_path('r-meta'))
                coloc_job.cpu(job_cpu)
                coloc_job.call(
                    coloc_runner,
                    hg38_map_chr_start_end,
                    eqtl_results_file,
                    celltype,
                    pheno_output_name,
                )
                manage_concurrency_for_job(coloc_job)

            else:
                print('No cis results for ' + gene + ' exist: skipping....')

    b.run(wait=False)


if __name__ == '__main__':
    main()
