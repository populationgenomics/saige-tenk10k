# Statistical colocalisation with SAIGE-QTL results

Using the [COLOC](https://chr1swallace.github.io/coloc/index.html) package.

## Contents

* [coloc runner](coloc_runner.py)
* [one more coloc runner for UKBB trait specifically](coloc_ukbb_runner.py)
* shell script to automate running on multiple phenotypes
  * [disease GWAS catalog traits](multi_pheno_coloc_runner.sh)
  * [UKBB blood traits](multi_pheno_ukbb_coloc_runner.sh)


## Combiner

I am using [Hope's script](https://github.com/populationgenomics/sv-workflows/blob/main/str/coloc/coloc_results_parser.py) to concatenate coloc results.

To run (from my local copy of `populationgenomics/sv-workflows/blob/main/str/coloc`):

```bash
analysis-runner --dataset "tenk10k" \
    --description "Combine coloc results" \
    --access-level "standard" \
    --output-dir "saige-qtl/tenk10k-genome-2-3-eur/output_files/241210/coloc-snp-only/sig_genes_only/" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-bioheart-main-analysis/saige-qtl/tenk10k-genome-2-3-eur/output_files/241210/coloc-snp-only/sig_genes_only \
    --celltypes=CD4_TCM \
    --phenos=ibd_liu2023
```

## Traits considered

### Autoimmune conditions

* inflammatory bowel disease (IBD)
* rheumatoid arthritis (RA)
* systemic lupus erythematosus (SLE)
* IgA nephropathy
* steroid sensitive nephrotic syndrome

### Cancer

* breast
* colorectal
* lung
* lymphoma
* prostate
* Non-Hodgkin’s lymphoma
* lymphocytic leukemia

### neuroinflammatory conditions

* Alzheimer’s disease
* Parkinson’s disease71

### Other

* COVID-19
* 44 UKBB traits from Margoliash et al. (Gymrek)
