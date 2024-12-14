#!/usr/bin/env python3
# pylint: disable=import-error,line-too-long

"""
Hail Batch workflow to perform association tests using SAIGE-QTL

- given all input files already generated
    - pheno + cov file (from get_anndata.py)
    - cis window file (from get_anndata.py)
    - genotype file (from get_genotype_vcf.py)
    - VRE genotypes (from get_genotype_vcf.py)
- builds saige commands (just add to str)
- run single-variant test using SAIGE-QTL (execute Rscript from command line)
- aggregate & summarise results


To run:

analysis-runner \
  --config saige_assoc_test.toml \
  --description "SAIGE-QTL CV association pipeline" \
  --memory='32G' \
  --storage='50G' \
  --dataset "tenk10k" \
  --access-level "full" \
  --output-dir "saige-qtl/tenk10k-genome-2-3-eur/output_files/241210" \
   python3 summarise_saige_assoc.py  --pheno-cov-files-path=gs://cpg-tenk10k-main/saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/pheno_cov_files \
      --cis-window-files-path=gs://cpg-tenk10k-main/saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/cis_window_files \
      --vre-files-prefix=gs://cpg-tenk10k-main/saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/genotypes/vds-tenk10k-genome-2-0
"""

import click
import json
import logging

from datetime import datetime
from os import getenv

from google.cloud import storage
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config, image_path, output_path, try_get_ar_guid
from cpg_utils.hail_batch import dataset_path, get_batch


def apply_job_settings(job: hb.batch.job.Job, job_name: str):
    """
    Apply settings to a job based on the name

    Args:
        job (hb.batch.job): job to apply settings to
        job_name (str): name used to find settings
    """
    # if this job has a section in config
    if job_settings := get_config()['saige']['job_specs'].get(job_name):
        # if memory setting in config - apply it
        if memory := job_settings.get('memory'):
            job.memory(memory)
        # if storage setting in config - apply it
        if storage := job_settings.get('storage'):
            job.storage(storage)
        # if cpu setting in config - apply it
        if cpu := job_settings.get('cpu'):
            job.cpu(cpu)


def summarise_cv_results(
    celltype: str,
    gene_results_path: str,
    summary_output_path: str,
):
    """
    Summarise gene-specific results
    """
    import logging
    import pandas as pd
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    existing_cv_assoc_results = [
        str(file)
        for file in to_path(gene_results_path).glob(f'*/{celltype}_*_cis_gene_pval')
    ]
    results_all_df = pd.concat(
        [
            pd.read_csv(to_path(pv_df), index_col=0)
            for pv_df in existing_cv_assoc_results
        ]
    )
    result_all_filename = to_path(output_path(summary_output_path, category='analysis'))
    logging.info(f'Write summary results to {result_all_filename}')
    with result_all_filename.open('w') as rf:
        results_all_df.to_csv(rf)


def create_second_job(vcf_path: str) -> hb.batch.job.Job:
    """
    Create a second job to run the single variant test
    """
    # get the size of the vcf file
    storage_client = storage.Client()
    bucket, filepath = vcf_path.removeprefix('gs://').split('/', 1)
    blob = storage_client.bucket(bucket).blob(filepath)
    blob.reload()  # refresh the blob to get the metadata
    size = blob.size // (1024**3)  # bytes to GB

    second_job = get_batch().new_job(name="saige-qtl part 2")
    apply_job_settings(second_job, 'sv_test')

    # VCF size, plus a 5GB buffer
    second_job.storage(f'{size + 5 }Gi')
    second_job.image(image_path('saige-qtl'))
    return second_job


@click.option(
    '--pheno-cov-files-path',
    default=dataset_path('saige-qtl/input_files/pheno_cov_files'),
)
@click.option(
    '--cis-window-files-path',
    default=dataset_path('saige-qtl/input_files/cis_window_files'),
)
@click.option(
    '--vre-files-prefix', default=dataset_path('saige-qtl/input_files/genotypes')
)
@click.option(
    '--writeout-file-prefix', default=dataset_path('saige-qtl', category='analysis')
)
@click.command()
def main(
    # outputs from gene_expression processing
    pheno_cov_files_path: str,
    cis_window_files_path: str,
    vre_files_prefix: str,
    # write out inputs and flags used for this run
    writeout_file_prefix: str,
):
    """
    Run SAIGE-QTL pipeline for all cell types
    """

    batch = get_batch('SAIGE-QTL pipeline')

    # define writeout file by type of pipeline and date and time
    date_and_time = datetime.today().strftime('%Y-%m-%d_%H:%M:%S')
    writeout_file = (
        f'{writeout_file_prefix}/saige_qtl_common_variant_pipeline_{date_and_time}.json'
    )

    # pull principal args from config
    celltypes: list[str] = get_config()['saige']['celltypes']
    celltype_jobs: dict[str, list] = dict()

    vre_plink_path = f'{vre_files_prefix}/vre_plink_2000_variants'

    # populate all the important params into a file for long-term reference
    writeout_dict: dict = {
        'ar_guid': try_get_ar_guid() or 'UNKNOWN',
        'results_output_path': output_path(''),
        'vre_plink_files_prefix_used': vre_plink_path,
        'pheno_cov_files_path_used': pheno_cov_files_path,
        'cis_window_files_path_used': cis_window_files_path,
        'saige_params': get_config()['saige'],
        'runtime_config': getenv('CPG_CONFIG_PATH'),
    }

    # summarise results (per cell type)
    for celltype in celltypes:
        logging.info(f'start summarising results for {celltype}')
        summary_output_path = (
            f'summary_stats/{celltype}_all_cis_cv_gene_level_results.tsv'
        )

        summarise_job = get_batch().new_python_job(
            f'Summarise CV results for {celltype}'
        )
        if celltype in celltype_jobs:
            summarise_job.depends_on(*celltype_jobs[celltype])
        summarise_job.call(
            summarise_cv_results,
            celltype=celltype,
            gene_results_path=output_path(celltype),
            summary_output_path=summary_output_path,
        )

    # write the file containing all important input parameters
    with to_path(writeout_file).open('wt') as out_handle:
        json.dump(writeout_dict, fp=out_handle, indent=4)

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
