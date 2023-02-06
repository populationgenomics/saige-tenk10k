# Hail batch workflow to run SAIGE on TenK10K data 

This is a hail batch pipeline to run the new [QTL version of SAIGE](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between rare genetic variants and single-cell gene expression from blood.

* **Plan A**: at present, just adapting our [CellRegMap Hail batch pipeline](https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py) hoping it can R / external code smoothly
* **Plan B**: if that fails, we may need to adapt [Konrad K's UKBB exomes analysis github](https://github.com/Nealelab/ukb_exomes), underlying [this paper](https://www.sciencedirect.com/science/article/pii/S2666979X22001100), or at least using [these python wrappers for SAIGE](https://github.com/Nealelab/ukb_common/blob/master/utils/saige_pipeline.py).

## Plan A

### preprocessing
* Hail query to filter WGS object to QC-passing, non ref-ref, rare (freq<5%) variants
* SAIGE R script to create sparse GRM
  * just once for all individuals, figure out which variants to use (just common? Ask Wei)

### gene-specific (and cell-type specific?)
* Hail query to filter object to relevant variants (gene-specific, maybe cell-type-specific also)
* Python (R?) to manipulate data to generate appropriate input files
* SAIGE running
  * null model (one per gene)
  * association (multiple?)

### cell-type-specific (all genes), aggregate
* back to python to aggregate results

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (1,000 individuals) - already generated
* Phase 1: OneK1K + BioHEART (2,000 individuals) - WGS done, scRNA-seq in progress
