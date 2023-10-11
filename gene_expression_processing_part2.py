#!/usr/bin/env python3


"""
PART 2 of 2 of the 'Gene expression processing workflow'

This script will:

- build chromosome & cell type specific phenotype covariate files
- create cis-window files for selected genes
- export pheno_cov files as tsv files
- export cis_window files as tsv files

More details in README
output files in tob_wgs_genetics/saige_qtl/input

    analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files" \
    --access-level "test" \
    --output-dir "tob_wgs_genetics/saige_qtl/hope-test-input" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
     gene_expression_processing_part2.py  --celltypes=B_IN --chromosomes=chr22 \
    --sample-mapping-file-path=gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv \
    --expression-files-prefix=hope-test

"""

import logging
import sys

import click
import hail as hl
import pandas as pd
from cpg_workflows.batch import get_batch
import hailtop.batch as hb
import scanpy
import json
import numpy as np
import hail as hl


from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path

from cpg_utils.hail_batch import dataset_path, output_path, init_batch

config = get_config()

SCANPY_IMAGE = config['images']['scanpy']

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr,
)

def get_celltype_covariates(
    expression_files_prefix: str,
    cell_type: str,
):
    """Obtain cell type specific covariates

    Input:
    - cell type of interest
    - covariate files prefix

    Output: covariate df for cell type of interest
    """
    # FIX the code to get cell type specific covariates

    #code to generate the covariate df file on line 91
    #expression_h5ad_path = to_path(
    #    dataset_path(
    #       f'scrna-seq/CellRegMap_input_files/expression_objects/sce22.h5ad'
    #    )
    # ).copy('here.h5ad')
    #adata = scanpy.read(expression_h5ad_path)
    #df = pd.DataFrame(adata.obsm['X_pca'], index = list(adata.obs.index)) #creates Pandas DF containing the PC loadings for each cell, index = CELL ID
    #df['OneK1K_ID'] = adata.obs['individual'] #adds in the Individual ID to the df
    #df.reset_index(inplace=True) #makes the CELL ID into a column
    #df_c = df[['index', 'OneK1K_ID', 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]] #select the first 15 PCs
    #df_c = df_c.rename(columns={0: "pce1", 1: "pce2",2: "pce3",3: "pce4",4: "pce5",5: "pce6",6: "pce7",7: "pce8",8: "pce9",9: "pce10",10: "pce11",11: "pce12",12: "pce13",13: "pce14",14: "pce15", "index":"Cell_ID"})

    #978 rows x 7 cols ##FYI to Anna - we lose some individuals here because 978 != number of individuals in df_c
    #covariates = pd.read_csv("gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/covariates_files/covariates.tsv", sep = "\t")

    #return pd.merge(df_c, covariates, left_on = "OneK1K_ID", right_on= "sampleid", how = "inner").drop('sampleid', axis=1) #we lose some individuals (as above), inner join - remove any NA fields

    covs_tsv_path = dataset_path('tob_wgs_genetics/saige_qtl/input/covariate_chr22_B_IN_tester_final.csv')
    covs_df = pd.read_csv(covs_tsv_path, sep=',', index_col=0)
    return covs_df


def build_pheno_cov_filename(
    gene_name, expression_files_prefix, celltype,chromosome, sample_mapping_file_path, ofile_path: str
):
    """
    Combine files to build final input
    reduce to one file per gene

    Input:
    - Expression anndata (filtered)
    - Sample mapping file path

    """
    #read in filtered anndata
    filtered_h5ad_path = to_path((output_path(f'filtered_{celltype}_{chromosome}.h5ad'))).copy('here.h5ad')
    expression_adata = scanpy.read(filtered_h5ad_path)

    #read in sample mapping file and covariate files
    smf_df = pd.read_csv(sample_mapping_file_path, sep='\t')
    cov_df = get_celltype_covariates(expression_files_prefix, celltype)

    ## TO FIX - we are not using raw data anymore

    #subset anndata by gene
    gene_index = np.where(expression_adata.raw.var.index == gene_name)[0][0]
    gene_mat = expression_adata.raw.X[:,gene_index].todense()

    expression_df = pd.DataFrame(
    data=gene_mat,
    index=expression_adata.obs.index,  #cell IDs
    columns=expression_adata.var_names[np.where(expression_adata.var_names == gene_name)] #complicated way to extract gene name
)
    expression_df['OneK1K_ID'] = expression_adata.obs['individual']

    # Make the index (cell IDs) into a column
    expression_df.reset_index(inplace=True)

    # Rename the columns
    expression_df.columns = ['Cell_ID', f'{gene_name}_raw_count', 'OneK1K_ID']

    #Merge data frame with sample mapping info and covariate file
    expression_df = pd.merge(expression_df, smf_df, on='OneK1K_ID', how='left')
    pheno_cov_df = pd.merge(expression_df, cov_df, on = "Cell_ID", how = "inner")

    #write to temp file path
    pheno_cov_df.to_csv(str(ofile_path), index=False)

def gene_info(x):
    # Extract ENSG and gene_level of evidence
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    g_id = list(filter(lambda x: 'gene_id' in x,  x.split(";")))[0].split("=")[1]
    g_id = g_id.split('.')[0] #removes the version number from ENSG ids
    return (g_name,g_id)

def get_gene_cis_file(chromosome:str, gene: str, window_size: int, ofile_path: str):
    """Get gene cis window file"""
    gencode = pd.read_table("gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/gencode.v42.annotation.gff3.gz", comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
    gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)
    gencode_genes["gene_name"],gencode_genes["ENSG"] = zip(*gencode_genes.attribute.apply(gene_info))

    #subset gencode annotation file for relevant chromosome
    gencode_genes = gencode_genes[gencode_genes['seqname']==chromosome]

    gene_info_gene = gencode_genes[gencode_genes['gene_name'] == gene]

    init_batch()
    # get chromosome, coordinates
    chrom = gene_info_gene['seqname'].values[0]
    start_coordinate = gene_info_gene['start'].values[0]
    end_coordinate = gene_info_gene['end'].values[0]
    # get gene body position (start and end) and add window
    left_boundary = max(1, start_coordinate - window_size)
    right_boundary = min(
        end_coordinate + window_size,
        hl.get_reference('GRCh38').lengths[chrom]
    )
    data = {'chromosome': chrom, 'start': left_boundary, 'end': right_boundary}
    pd.DataFrame(data, index=[gene]).to_csv(str(ofile_path))


@click.command()
@click.option('--celltypes')
@click.option('--chromosomes')
@click.option('--expression-files-prefix')
@click.option('--sample-mapping-file-path')
@click.option('--min-pct-expr', type=int, default =5)
@click.option('--cis-window-size', type=int,default=100000)
@click.option(
    '--max-gene-concurrency',
    type=int,
    default=50,
    help=(
        'To avoid resource starvation, set this concurrency to limit '
        'horizontal scale. Higher numbers have a better walltime, but '
        'risk jobs that are stuck (which are expensive)'
    ),
)
def main(
    celltypes: str,
    chromosomes: str,
    expression_files_prefix: str,
    sample_mapping_file_path: str,
    cis_window_size: int,
    max_gene_concurrency=int
):
    """
    Run expression processing pipeline
    """
    config = get_config()
    b = get_batch()
    init_batch()


    logging.info(f'Cell types to run: {celltypes}')
    logging.info(f'Chromosomes to run: {chromosomes}')

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_gene_concurrency:
            job.depends_on(_dependent_jobs[-max_gene_concurrency])
        _dependent_jobs.append(job)

    for celltype in celltypes.split(','):
        for chromosome in chromosomes.split(','):

            #load array containing the genes from filtered anndata
            with to_path(output_path(f'{chromosome}_{celltype}_filtered_genes.json')).open('r') as read_handle:
                genes = json.load(read_handle)

            # combine files for each gene
            # pylint: disable=no-member
            for gene in genes:
                pheno_cov_job = b.new_python_job(name=f'Build phenotype-covariate files for {gene} [{celltype};{chromosome}]')
                pheno_cov_job.storage('8G')
                pheno_cov_job.cpu(2)
                pheno_cov_job.image(config['workflow']['driver_image'])
                pheno_cov_job.call(
                    build_pheno_cov_filename,gene,"", celltype, chromosome, sample_mapping_file_path,pheno_cov_job.ofile)
                pheno_cov_job.ofile.add_extension('.csv')
                b.write_output(pheno_cov_job.ofile, output_path(
                      f'input_files/pheno_cov_files/{gene}_{celltype}.csv'
                    ))
                manage_concurrency_for_job(pheno_cov_job)

                # add gene cis window file
                gene_cis_job = b.new_python_job(name=f'Build cis window files for {gene} [{celltype};{chromosome}]')
                gene_cis_job.image(config['workflow']['driver_image'])
                gene_cis_job.call(
                    get_gene_cis_file,chromosome,gene,cis_window_size,gene_cis_job.ofile
                )
                b.write_output(gene_cis_job.ofile, output_path(f'input_files/cis_window_files/{gene}_{cis_window_size}bp.csv'))

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

