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
    --celltypes "B_naive"

"""

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, image_path, output_path


def coloc_runner(gwas, eqtl_file_path, celltype, coloc_results_file):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

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
    eqtl = pd.read_csv(
        eqtl_file_path,
        sep='\t',
    )
    eqtl['beta'] = eqtl['BETA']
    eqtl['se'] = eqtl['SE']
    eqtl['position'] = eqtl['POS']
    eqtl['snp'] = eqtl['MarkerID'].str.replace(':', '_', regex=False)
    # while I figure out if it's easy to extract sdY, give N and MAF instead
    # https://chr1swallace.github.io/coloc/articles/a02_data.html#what-if-i-dont-have-sdy
    eqtl['MAF'] = eqtl['AF_Allele2'].apply(lambda af: min(af, (1 - af)))
    gene = eqtl_file_path.split('/')[-1].replace(f'{celltype}_', '').replace('_cis', '')
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
    eqtl_r = eqtl_r %>% select(beta, varbeta, position, snp, N, MAF)

    eqtl_r = eqtl_r %>% as.list()
    eqtl_r$type = 'quant'

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
        coloc_results_file,
        sep='\t',
        index=False,
    )


@click.option(
    '--egenes-files-path',
    help='Path to the gene-level summary files',
    default='gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files/summary_stats',
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
@click.option(
    '--celltypes', help='Cell types to run, single str, comma separated', default='ASDC'
)
@click.option('--cis-window-size', help='Cis window size used', default=100000)
@click.option('--fdr-threshold', help='FDR threshold', default=0.05)
@click.option(
    '--max-parallel-jobs', help='Maximum number of parallel jobs to run', default=500
)
@click.option(
    '--pheno-output-name', help='Phenotype output name', default='covid_GCST011071'
)
@click.option('--job-cpu', help='Number of CPUs to use for each job', default=0.25)
@click.command()
def main(
    snp_cis_dir: str,
    egenes_files_path: str,
    celltypes: str,
    snp_gwas_file: str,
    pheno_output_name: str,
    max_parallel_jobs: int,
    cis_window_size: int,
    fdr_threshold: float,
    job_cpu: float,
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
        'gs://cpg-bioheart-test/saige-qtl/300-libraries/combined_anndata_obs_vars/300_libraries_concatenated_harmony_filtered_vars.csv',
    )
    hg38_map = pd.read_csv(
        snp_gwas_file,
        sep='\t',
    )

    b = get_batch(name=f'Run coloc:{pheno_output_name}')

    for celltype in celltypes.split(','):

        # read in eGenes file
        egenes_file = (
            egenes_file
        ) = f'{egenes_files_path}/{celltype}_all_cis_cv_gene_level_results.tsv'
        result_df_cfm = pd.read_csv(
            egenes_file,
            sep='\t',
        )
        result_df_cfm['chr'] = result_df_cfm['top_MarkerID'].apply(
            lambda snp: 'chr' + snp.split(':')[0]
        )
        result_df_cfm = result_df_cfm[
            result_df_cfm['ACAT_p'] < fdr_threshold
        ]  # filter for sc-eQTLs with p-value < fdr_threshold

        for gene in result_df_cfm['gene']:
            chrom = result_df_cfm[result_df_cfm['gene'] == gene]['chr'].iloc[0]
            coloc_results_file = output_path(
                f"coloc-snp-only/sig_genes_only/{pheno_output_name}/{celltype}/{gene}_{cis_window_size}.tsv",
                'analysis',
            )
            if to_path(coloc_results_file).exists():
                continue

            eqtl_results_file = (
                f'{snp_cis_dir}/{celltype}/{chrom}/{celltype}_{gene}_cis'
            )
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
                if hg38_map_chr_start_end['p_value'].min() > 5e-8:
                    print(
                        'No significant SNP GWAS data for '
                        + gene
                        + ' in the cis-window: skipping....'
                    )
                    continue
                print('Extracted SNP GWAS data for ' + gene)

                # run coloc
                coloc_job = b.new_python_job(
                    f'Coloc for {gene}: {celltype}',
                )
                coloc_job.image(image_path('r-meta'))
                coloc_job.cpu(job_cpu)
                coloc_job.call(
                    coloc_runner,
                    hg38_map_chr_start_end,
                    eqtl_results_file,
                    celltype,
                    coloc_results_file,
                )
                manage_concurrency_for_job(coloc_job)

            else:
                print('No cis results for ' + gene + ' exist: skipping....')

    b.run(wait=False)


if __name__ == '__main__':
    main()
