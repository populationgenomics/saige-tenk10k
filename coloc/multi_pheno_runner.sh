#!/bin/bash

# Define the cell types as a comma-separated string
celltypes='B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,gdT,dnT,CD4_TCM'
# celltypes='B_naive' # test

# Convert the cell types into an array
IFS=',' read -ra celltype_array <<< "$celltypes"

# Loop through each phenotype
for i in \
    "alzheimer_GCST90027158 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90027158.h_parsed.tsv" \
    "breastca_GCST004988 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004988.h_parsed.tsv" \
    "colorectalca_GCST90129505 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/colorectalca_GCST90129505_parsed.tsv" \
    "covid_GCST011071 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv" \
    "ibd_liu2023 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/ibd_EAS_EUR_SiKJEF_meta_IBD.tsv" \
    "NHL_GCST90011819 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/NHL_GCST90011819_parsed.tsv" \
    "lungca_GCST004748 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004748.h_parsed.tsv" \
    "lymphoma_GCST90018878 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90018878.h_parsed.tsv" \
    "parkinson_GCST009325 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST009325.h_parsed.tsv" \
    "prostateca_GCST90274713 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90274713.h_parsed.tsv" \
    "ra_GCST90132223 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90132223_parsed.tsv"\
    "sle_GCST003156 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/bentham_2015_26502338_sle_parsed.tsv'"\
    "myeloproliferative_GCST90000032 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/myeloproliferative_GCST90000032_parsed.tsv" \
    "lymphocytic_leukemia_GCST90011814 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/lymphocytic_leukemia_GCST90011814_parsed.tsv" \
    "nephrotic_GCST90258619 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/nephrotic_GCST90258619_parsed.tsv" \
    "kiryluk_IgAN gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/Kiryluk_IgAN_parsed.tsv" \
    "trujillo_methylation_eQTLs gs://cpg-bioheart-test/str/Trujillo_methylation_eQTLs/hg38_STRs_SNVs_parsed.tsv" \
    ; do
    set -- $i
    # Loop through each cell type
    for celltype in "${celltype_array[@]}"; do
        echo "Running coloc analysis for cell type: $celltype, phenotype: $1"
        analysis-runner --dataset "bioheart" \
        --description "Run coloc for eGenes identified by SAIGE-QTL analysis" \
        --access-level "full" \
        --memory='8G' \
        --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
        --output-dir "saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/" \
        coloc/coloc_runner.py \
        --egenes-files-path=gs://cpg-bioheart-main-analysis/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/output_files/summary_stats \
        --snp-cis-dir=gs://cpg-bioheart-main-analysis/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/output_files \
        --snp-gwas-file="$2" \
        --pheno-output-name="$1" \
        --celltypes="$celltype"
    done
done


## test

# # Loop through each phenotype
# for i in \
#     "alzheimer_GCST90027158 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90027158.h_parsed.tsv" \
#     "breastca_GCST004988 gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004988.h_parsed.tsv" \
#     ; do
#     set -- $i
#     # Loop through each cell type
#     for celltype in "${celltype_array[@]}"; do
#         echo "Running coloc analysis for cell type: $celltype, phenotype: $1"
#         analysis-runner --dataset "bioheart" \
#         --description "Run coloc for eGenes identified by SAIGE-QTL analysis" \
#         --access-level "test" \
#         --memory='8G' \
#         --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
#         --output-dir "saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/" \
#         coloc/coloc_runner.py \
#         --egenes-files-path=gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files/summary_stats \
#         --snp-cis-dir=gs://cpg-bioheart-test-analysis/saige-qtl/bioheart_n990_and_tob_n1055/241004_n100/output_files \
#         --snp-gwas-file="$2" \
#         --pheno-output-name="$1" \
#         --celltypes="$celltype" --fdr-threshold=1
#     done
# done
