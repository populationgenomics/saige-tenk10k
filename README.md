# Hail batch workflow to run SAIGE-QTL on TenK10K data

This is a hail batch pipeline to run the new [QTL version of SAIGE](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between common and rare genetic variants and single-cell gene expression from blood.
First, this will be run on the TOB (AKA OneK1K) and then BioHEART datasets as part of the TenK10K project (see *Data* below), but in the future all datasets within OurDNA (with scRNA-seq + WGS data) will go through this pipeline as well.

The pipeline is split into four parts, to make for more flexible usage:

* Genotype processing (SNV): this involves sample and variant QC of the WGS data, and genotype file preparation specifically for common and rare single-nucleotide variants
* Expression (phenotype) processing: this involves processing of the scRNA-seq data, and preparation of the pheno_cov file
* Input files check and preparation: this involves combining data from above and cross-check for consistency prior to running SAIGE-QTL
* Association testing: prepare and run SAIGE-QTL commands for association mapping


## Genotypes preprocessing (once per cohort, e.g. TOB)

Function name: filter_variants.

Hail query to filter WGS object to

* samples that are: i) QC-passing, ii) present in the scRNA-seq dataset
* variants that re: i) QC-passing, ii) non ref-ref variants, and iii) (for now) indels and multi-allelic SNPs.

It outputs three objects:

* MT object, all retained samples and variants
* plink object for only 2,000 variants (MAC>20), after LD pruning - this is for the estimation of the variance ratio (VR plinks)

consider ouputting one single MT object instead.

### To run

```bash
analysis-runner \
    --dataset tob-wgs \
    --access-level standard \
    --output-dir 'tob_wgs_genetics/saige_qtl/input' \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:dev \
    --description 'WGS processing batch job' \
    python3 genotypes_processing.py \
      --sample-mapping-file-tsv 'tob_wgs_genetics/saige_qtl/input/smf.tsv'
```

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (1,000 individuals) - already generated (both WGS and scRNA-seq, though with an older technology)
* Phase 1: OneK1K + BioHEART (2,000 individuals) - WGS done, scRNA-seq in progress (new kit, which should result in many more cells per individual)
* Phase 2/final: aim is ~ 10,000 individuals from the (extended) TOB/OneK1K, BioHEART, ADAPT, LBIO and AIM cohorts (nothing generated besides current stage of Phase1)
