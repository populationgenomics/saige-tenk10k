#!/bin/bash

# Define the cell types as a comma-separated string
celltypes='B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,gdT,dnT,CD4_TCM'

# Convert the cell types into an array
IFS=',' read -ra celltype_array <<< "$celltypes"

# Define phenotypes and file paths as parallel arrays
pheno_names=(
    "alzheimer_GCST90027158"
    "breastca_GCST004988"
    "colorectalca_GCST90129505"
    "covid_GCST011071"
    "ibd_liu2023"
    "NHL_GCST90011819"
    "lungca_GCST004748"
    "lymphoma_GCST90018878"
    "parkinson_GCST009325"
    "prostateca_GCST90274713"
    "ra_GCST90132223"
    "sle_GCST003156"
    "myeloproliferative_GCST90000032"
    "lymphocytic_leukemia_GCST90011814"
    "nephrotic_GCST90258619"
    "kiryluk_IgAN"
)

pheno_files=(
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90027158.h_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004988.h_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/colorectalca_GCST90129505_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/ibd_EAS_EUR_SiKJEF_meta_IBD.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/NHL_GCST90011819_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004748.h_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90018878.h_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST009325.h_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90274713.h_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90132223_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/bentham_2015_26502338_sle_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/myeloproliferative_GCST90000032_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/lymphocytic_leukemia_GCST90011814_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/nephrotic_GCST90258619_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/Kiryluk_IgAN_parsed.tsv"
)

# Loop through each cell type and phenotype
for celltype in "${celltype_array[@]}"; do
    for i in "${!pheno_names[@]}"; do
        pheno="${pheno_names[$i]}"
        filepath="${pheno_files[$i]}"
        analysis-runner --dataset "bioheart" \
        --description "Run coloc for eGenes identified by STR analysis" \
        --access-level "test" \
        --memory='16G' \
        --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
        --output-dir "tenk10k/str/associatr/final_freeze" \
        coloc_runner.py \
        --snp-gwas-file="$filepath" \
        --pheno-output-name="$pheno" \
        --celltypes="$celltype" \
        --max-parallel-jobs 10000
    done
done
