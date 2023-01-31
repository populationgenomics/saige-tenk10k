# Hail batch workflow to run SAIGE on TenK10K data 

This is a hail batch pipeline to run the new [QTL version of SAIGE](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between rare genetic variants and single-cell gene expression from blood.

* **Plan A**: at present, just adapting our [CellRegMap Hail batch pipeline](https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py) hoping it can R / external code smoothly
* **Plan B**: if that fails, we may need to adapt [Konrad K's UKBB exomes analysis github](https://github.com/Nealelab/ukb_exomes), underlying [this paper](https://www.sciencedirect.com/science/article/pii/S2666979X22001100), or at least using [these python wrappers for SAIGE](https://github.com/Nealelab/ukb_common/blob/master/utils/saige_pipeline.py).

## Plan A

* Hail query to filter WGS object to QC-passing, non ref-ref, rare (freq<5%) variants
* Hail query to filter object to relevant variants (gene-specific, maybe cell-type-specific also)
* Python (R?) to manipulate data to generate appropriate input files
* SAIGE running
  * create SPARSE GRM?
  * null model (one per gene)
  * association (multiple?)
* back to python to aggregate results
