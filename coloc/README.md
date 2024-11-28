# Statistical colocalisation with SAIGE-QTL results

Using COLOC {ref}

## Contents

* coloc runner
* maybe one more coloc runner for UKBB traits
* shell script to automate running on multiple phenotypes

## Combiner

I am using [Hope's script](https://github.com/populationgenomics/sv-workflows/blob/main/str/coloc/coloc_results_parser.py) to concatenate coloc results.

To run (from my local copy of `populationgenomics/sv-workflows/blob/main/str/coloc`):

```bash
analysis-runner --dataset "bioheart" \
    --description "Combine coloc results" \
    --access-level "standard" \
    --output-dir "saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/coloc-snp-only/sig_genes_only/" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-bioheart-main-analysis/saige-qtl/bioheart_n787_and_tob_n960/241008_ashg/coloc-snp-only/sig_genes_only \
    --celltypes=CD4_TCM \
    --phenos=ibd_liu2023
```

## Traits considered

alzheimer_GCST90027158 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90027158.h_parsed.tsv'
breastca_GCST004988 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004988.h_parsed.tsv'
colorectalca_GCST90129505 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/colorectalca_GCST90129505_parsed.tsv'
covid_GCST011071 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv'
ibd_liu2023 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/ibd_EAS_EUR_SiKJEF_meta_IBD.tsv'
NHL_GCST90011819 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/NHL_GCST90011819_parsed.tsv'
lungca_GCST004748 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004748.h_parsed.tsv'
lymphoma_GCST90018878 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90018878.h_parsed.tsv'
parkinson_GCST009325 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST009325.h_parsed.tsv'
prostateca_GCST90274713 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90274713.h_parsed.tsv'
ra_GCST90132223 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90132223_parsed.tsv'
sle_GCST003156 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/bentham_2015_26502338_sle_parsed.tsv'
myeloproliferative_GCST90000032 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/myeloproliferative_GCST90000032_parsed.tsv'
lymphocytic_leukemia_GCST90011814 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/lymphocytic_leukemia_GCST90011814_parsed.tsv'
nephrotic_GCST90258619 = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/nephrotic_GCST90258619_parsed.tsv'
kiryluk_IgAN = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/Kiryluk_IgAN_parsed.tsv'
trujillo_methylation_eQTLs = 'gs://cpg-bioheart-test/str/Trujillo_methylation_eQTLs/hg38_STRs_SNVs_parsed.tsv'

44 UKBB traits from Margoliash et al. (Gymrek)
'gymrek-ukbb-alanine-aminotransferase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_alanine_aminotransferase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-albumin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_albumin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-alkaline_phosphatase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_alkaline_phosphatase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-apolipoprotein_a': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_apolipoprotein_a_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-apolipoprotein_b': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_apolipoprotein_b_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-aspartate_aminotransferase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_aspartate_aminotransferase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-c_reactive_protein': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_c_reactive_protein_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-calcium': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_calcium_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-cholesterol': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_cholesterol_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-creatinine': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_creatinine_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-cystatin_c': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_cystatin_c_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-eosinophil_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_eosinophil_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-eosinophil_percent': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_eosinophil_percent_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-gamma_glutamyltransferase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_gamma_glutamyltransferase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-glucose': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_glucose_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-glycated_haemoglobin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_glycated_haemoglobin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-haematocrit': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_haematocrit_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-haemoglobin_concentration': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_haemoglobin_concentration_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-hdl_cholesterol': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_hdl_cholesterol_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-igf_1': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_igf_1_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-ldl_cholesterol_direct': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_ldl_cholesterol_direct_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-phosphate': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_phosphate_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-shbg': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_shbg_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-total_bilirubin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_total_bilirubin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-total_protein': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_total_protein_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-triglycerides': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_triglycerides_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-urate': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_urate_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-urea': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_urea_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-vitamin_d': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_vitamin_d_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-white_blood_cell_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_white_blood_cell_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-lymphocyte_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_lymphocyte_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-lymphocyte_percent': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_lymphocyte_percent_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_corpuscular_haemoglobin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_corpuscular_haemoglobin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_corpuscular_haemoglobin_concentration': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_corpuscular_haemoglobin_concentration_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_corpuscular_volume': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_corpuscular_volume_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_platelet_volume': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_platelet_volume_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_sphered_cell_volume': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_sphered_cell_volume_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-neutrophil_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_neutrophil_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-neutrophil_percent': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_neutrophil_percent_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-platelet_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_platelet_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-platelet_crit': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_platelet_crit_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-platelet_distribution_width': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_platelet_distribution_width_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-red_blood_cell_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_red_blood_cell_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-red_blood_cell_distribution_width': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_red_blood_cell_distribution_width_snp_str_gwas_results_hg38_{chrom}.tab.gz',


