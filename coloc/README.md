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