#!/usr/bin/env python3

"""
Copy of https://github.com/populationgenomics/sv-workflows/blob/main/str/coloc/coloc_ukbb_runner.py
to rework for common variant SAIGE-QTL results


This script performs SNP+STR colocalisation analysis betweeen eGenes identified by pseudobulk STR analysis and UKBB GWAS signals.
Assumes that the SNP GWAS data has been pre-processed with the following columns: 'chromosome', 'position' (hg38 bp), 'snp'(chromosome_position_refallele_effectallele), 'beta', 'varbeta'
1) Extract eGenes (FDR 5%) and those where STR pval < 5e-8
2) Extract the SNP GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene (if the GWAS data has at least one variant with pval <5e-8)
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by STR analysis" \
    --access-level "test" \
    --memory='4G' \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr" \
    coloc_ukbb_runner.py \
    --pheno-output-name=gymrek-ukbb-apolipoprotein-a \
    --celltypes "CD16_Mono" \
    --max-parallel-jobs 10000

"""

import gzip

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
    gwas_r$type = 'quant'
    gwas_r$sdY = 1

    ''',
    )
    eqtl = pd.read_csv(
        eqtl_file_path,
        sep='\t',
    )
    gene = eqtl_file_path.split('/')[-1].split('_')[0]
    eqtl['beta'] = eqtl['coeff_meta']
    eqtl['se'] = eqtl['se_meta']
    eqtl['position'] = eqtl['pos']
    eqtl['snp'] = eqtl['chr'] + '_' + eqtl['position'].astype(str) + '_' + eqtl['motif']
    eqtl['snp'] = eqtl['snp'].str.replace('-', '_', regex=False)

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
        output_path(f"coloc/sig_str_and_gwas_hit/{pheno_output_name}/{celltype}/{gene}_100kb.tsv", 'analysis'),
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
    default='gs://cpg-bioheart-test/str/associatr/snps_and_strs/tob_n1055_and_bioheart_n990/meta_results/meta_results',
)
@click.option('--celltypes', help='Cell types to run', default='ASDC')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs to run', default=500)
@click.option('--pheno-output-name', help='Phenotype output name', default='covid_GCST011071')
@click.option('--job-cpu', help='Number of CPUs to use for each job', default=0.25)
@click.command()
def main(snp_cis_dir, egenes_file, celltypes, pheno_output_name, max_parallel_jobs, job_cpu):
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
    # read in eGenes file
    egenes = pd.read_csv(
        egenes_file,
        sep='\t',
        usecols=['chr', 'pos', 'pval_meta', 'motif', 'susie_pip', 'gene', 'finemap_prob', 'celltype', 'ref_len'],
    )

    result_df_cfm = egenes
    result_df_cfm['variant_type'] = result_df_cfm['motif'].str.contains('-').map({True: 'SNV', False: 'STR'})
    result_df_cfm_str = result_df_cfm[result_df_cfm['variant_type'] == 'STR']  # filter for STRs
    result_df_cfm_str = result_df_cfm_str[result_df_cfm_str['pval_meta'] < 5e-8]  # filter for STRs with p-value < 5e-8
    result_df_cfm_str = result_df_cfm_str.drop_duplicates(
        subset=['gene', 'celltype'],
    )  # drop duplicates (ie pull out the distinct genes in each celltype)
    result_df_cfm_str['gene'] = result_df_cfm_str['gene'].str.replace(
        '.tsv',
        '',
        regex=False,
    )  # remove .tsv from gene names (artefact of the data file)
    b = get_batch(name=f'Run coloc:{pheno_output_name} and {celltypes}')

    for celltype in celltypes.split(','):
        result_df_cfm_str_celltype = result_df_cfm_str[
            result_df_cfm_str['celltype'] == celltype
        ]  # filter for the celltype of interest
        for chrom in result_df_cfm_str_celltype['chr'].unique():
            result_df_cfm_str_celltype_chrom = result_df_cfm_str_celltype[result_df_cfm_str_celltype['chr'] == chrom]
            phenotype = pheno_output_name.split('-')[-1]
            chr_gwas_file = f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_{phenotype}_snp_str_gwas_results_hg38_{chrom}.tab.gz'

            with gzip.open(to_path(chr_gwas_file), 'rb') as f:
                hg38_map = pd.read_csv(f, sep='\t')

            for gene in result_df_cfm_str_celltype_chrom['gene']:
                if to_path(
                    output_path(
                        f"coloc/sig_str_and_gwas_hit/{pheno_output_name}/{celltype}/{gene}_100kb.tsv",
                        'analysis',
                    ),
                ).exists():
                    continue
                if to_path(f'{snp_cis_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv').exists():
                    print('Cis results for ' + gene + ' exist: proceed with coloc')

                    # extract the coordinates for the cis-window (gene +/- 100kB)
                    gene_table = var_table[var_table['gene_ids'] == gene]
                    start = float(gene_table['start'].astype(float)) - 100000
                    end = float(gene_table['end'].astype(float)) + 100000
                    hg38_map_chr = hg38_map[hg38_map['chromosome'] == (chrom)]
                    hg38_map_chr_start = hg38_map_chr[hg38_map_chr['position'] >= start]
                    hg38_map_chr_start_end = hg38_map_chr_start[hg38_map_chr_start['position'] <= end]
                    if hg38_map_chr_start_end.empty:
                        print('No GWAS data for ' + gene + ' in the cis-window: skipping....')
                        continue
                    # check if the p-value column contains at least one value which is <=5e-8:
                    if hg38_map_chr_start_end['p_value'].min() > 5e-8:
                        print('No significant SNP STR GWAS data for ' + gene + ' in the cis-window: skipping....')
                        continue
                    print('Extracted GWAS data for ' + gene)

                    # run coloc
                    coloc_job = b.new_python_job(
                        f'Coloc for {gene}: {celltype}: {phenotype}',
                    )
                    f'{snp_cis_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv'
                    coloc_job.image(image_path('r-meta'))
                    coloc_job.cpu(job_cpu)
                    coloc_job.call(
                        coloc_runner,
                        hg38_map_chr_start_end,
                        f'{snp_cis_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv',
                        celltype,
                        pheno_output_name,
                    )
                    manage_concurrency_for_job(coloc_job)

                else:
                    print('No cis results for ' + gene + ' exist: skipping....')

    b.run(wait=False)


if __name__ == '__main__':
    main()