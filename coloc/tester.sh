#!/bin/bash

# Define the cell types as a comma-separated string
celltypes='B_naive'

# Convert the cell types into an array
IFS=',' read -ra celltype_array <<< "$celltypes"

# Define phenotypes and file paths as parallel arrays
pheno_names=(
    "covid_GCST011071"
    "NHL_GCST90011819"
)

pheno_files=(
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv"
    "gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/NHL_GCST90011819_parsed.tsv"
)

# Loop through each cell type and phenotype
for celltype in "${celltype_array[@]}"; do
    for i in "${!pheno_names[@]}"; do
        pheno="${pheno_names[$i]}"
        filepath="${pheno_files[$i]}"
        analysis-runner --dataset "tenk10k" \
        --description "Run coloc for eGenes identified by SAIGE-QTL" \
        --access-level "test" \
        --memory='32G' \
        --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
        --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/" \
        coloc/coloc_runner.py \
        --egenes-files-path=gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files/summary_stats \
        --snp-cis-dir=gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files \
        --snp-gwas-file="$filepath" \
        --pheno-output-name="$pheno" \
        --celltypes="$celltype" \
        --max-parallel-jobs 10000 --fdr-threshold=1 --gwas-significance-threshold=0.05
    done
done