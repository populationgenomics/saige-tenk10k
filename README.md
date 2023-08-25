# Hail batch workflow to run SAIGE-QTL on TenK10K data

This is a hail batch pipeline to run the new [QTL version of SAIGE](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between both common and rare genetic variants and single-cell gene expression from blood.
First, this will be run on the TOB (AKA OneK1K) and then BioHEART datasets as part of phase 1 of the TenK10K project (see *Data* below), but in the future all datasets within OurDNA (with scRNA-seq + WGS data) will go through this pipeline as well.

The pipeline is split into four parts, to make for more flexible usage:

1. Genotype processing (SNV): this involves sample and variant QC of the WGS data, and genotype file preparation specifically for common and rare single-nucleotide variants
2. Expression (phenotype) processing: this involves processing of the scRNA-seq data, and preparation of the pheno_cov file
3. Input files check and preparation: this involves combining data from above and cross-check for consistency prior to running SAIGE-QTL
4. Association testing: prepare and run SAIGE-QTL commands for association mapping

Only [1] is ready for now.
Working on [2] in this branch.

## Genotypes preprocessing (once per cohort, e.g., TOB)

Function name: filter_variants.

Hail query to filter WGS object to

* samples that are: i) QC-passing, ii) present in the scRNA-seq dataset
* variants that are: i) QC-passing, ii) non ref-ref variants, and iii) (for now) indels and multi-allelic SNPs.

It outputs two objects:

* MT object, all retained samples and variants (both common & rare at this stage)
* plink object for only 2,000 variants (minor allele count>20), after LD pruning - this is for the estimation of the variance ratio (VR plinks)


### To run

```bash
analysis-runner \
    --dataset tob-wgs \
    --access-level standard \
    --output-dir 'tob_wgs_genetics/saige_qtl/input' \
    --description 'WGS processing batch job' \
    python3 genotypes_processing.py \
      --mt-path 'mt/v7.mt' \
      --sample-mapping-file-tsv 'scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
```

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (1,000 individuals) - already generated (both WGS and scRNA-seq, though with an older technology)
* Phase 1: OneK1K + BioHEART (2,000 individuals) - WGS done, scRNA-seq in progress (new kit, which should result in many more cells per individual)
* Phase 2/final: aim is ~ 10,000 individuals from the (extended) TOB/OneK1K, BioHEART, ADAPT, LBIO and AIM cohorts (nothing generated besides current stage of Phase 1)

## Additional resources

* [SAIGE-QTL pipeline flowchar GSlides](https://docs.google.com/presentation/d/1OhNiA6DaP9ZGlAbh8uZuZWzvrrr_QwvJwJ_lVPBZoic/edit#slide=id.g25daf727307_0_102)
* [SAIGE-QTL pipeline notes GDoc](https://docs.google.com/document/d/1t11VafeU1THA4X58keHd5oPVglTYiY3DKC7P05GHCzw/edit)
