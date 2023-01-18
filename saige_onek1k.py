#!/usr/bin/env python3
# pylint: disable=import-error

__author__ = 'annacuomo'

"""
Hail Batch workflow for the rare-variant association analysis, including:
- get relevant variants around a gene and export genotypes as plink files
- generate input files for association tests
- run association tests
"""

# python imports
import click
import logging

import hail as hl
import hailtop.batch as hb

# from cpg_utils import to_path
from cpg_utils.hail_batch import (
#     copy_common_env,
#     dataset_path,
    get_config,
#     init_batch,
#     output_path,
#     remote_tmpdir,
)


# import self scripts
from saige_onek1k import *

# Docker images
HAIL_DOCKER_IMAGE = ''  # would be good to use a fixed one (that doesn't always update, for reproducibility)
SAIGE_DOCKER_IMAGE = 'wzhou88/saige:1.1.6.3'  # from https://saigegit.github.io/SAIGE-doc/docs/Installation_docker.html

config = get_config()

# main function
@click.command()
@click.option('--celltypes')
def saige_pipeline(
    celltypes: str,
):
  sb = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
  batch = hb.Batch('CellRegMap pipeline', backend=sb)
  # set jobs running
  batch.run(wait=False)
  
  
if __name__ == '__main__':
    saige_pipeline()
