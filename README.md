# Hail batch workflow to run SAIGE on TenK10K data

This is a hail batch pipeline to run the new [QTL version of SAIGE](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between common and rare genetic variants and single-cell gene expression from blood.

* **Plan A**: at present, just adapting our [CellRegMap Hail batch pipeline](https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py) hoping it can run R / external code smoothly.
* **Plan B**: if that fails, we may need to adapt [Konrad K's UKBB exomes analysis github](https://github.com/Nealelab/ukb_exomes), underlying [this paper](https://www.sciencedirect.com/science/article/pii/S2666979X22001100), or at least using [these python wrappers for SAIGE](https://github.com/Nealelab/ukb_common/blob/master/utils/saige_pipeline.py).

## Plan A

### Genotypes preprocessing (once per cohort)

Hail query to filter WGS object to i) QC-passing, ii) non ref-ref variants, and considering only samples with scRNA-seq data.

It outputs three objects:

* MT object, rare (freq<5%) variants
* MT object, common (freq>1%) variants
* plink object for only 2,000 variants (MAC>20), after LD pruning - this is for the estimation of the variance ratio (VR plinks)

<!-- # skip for now - unrelated individuals
* SAIGE R script to create sparse GRM
  * just once for all individuals, all variants after LD-pruning, and MAF>1% -->

### Expression preprocessing (once per cell type)

Python script to combine expression (pheno), covariates into a single pheno_cov file as input to saige-qtl.

Inputs:

* chromosome-specific scanpy (AnnData) objects, single-cell expression for all cells, all genes for that chromosome (sctransformed sc-counts)
* covariates: combination of 1) exploded donor-level covariates (e.g. sex) and 2) cell-level covariates (PCs, batch)
* sample mapping file: matching donor ID (onek1k vs CPG) and cell barcodes, including cell type labels

Output file (one per cell type):

* text file concatenating covs, expression of genes, and individual id (same as plink files) for all cells from that cell type

### Variant selection (once per gene)

Hail query to filter object to relevant variants

* for common variants, only genomic proximity (+/-100kb)
* for rare variants, i) genomic proximity (+/-50kb), 2) regulatory consequences (vep), 3) open chromatin (any cell type)

Outputs:

* plink files (.bed, .bim, .fam) for common variants
* plink files (.bed, .bim, .fam) for rare variants
* group file for rare variants. For each region(=gene?)
  * row1: region1 | var    | 1:100010:A:C | 1:100016:A:C
  * row2: region1 | anno   | promoter     | enhancer
  * row3: region1 | weight | 0.5          | 0.3

### Run association (for each gene, cell type combination)

For both step 1 and step2:

  * Python job to build command ("Rscript --arg1 par1 --arg2 par2..)
  * Job to run command

#### Step1: Fit null model

Inputs:

* Pheno cov file
* VR plinks
* info on what gene, what covs, what likelihood to use

Outputs:

* model object (.rda)
* variance ratio estimate (.txt)

#### Step2: Run association

Inputs:

* model object (.rda)
* variance ratio estimate (.txt)
* genotypes to test (plink)
* rare only: group file

Output:

* txt file with all results

### Results aggregation (once per cell type)

Python to aggregate all results from Step2 above

### To run (this may need to be updated)

```bash
analysis-runner \
    --dataset tob-wgs \
    --access-level standard \
    --output-dir "tob_wgs_genetics/saige_qtl/" \
    --image ? \
    --description "Saige QTL batch job" \
    python3 saige_tenk10k.py \
      --input-files-prefix tob_wgs_genetics/saige_qtl/input \
      --sample-mapping-file-tsv tob_wgs_genetics/saige_qtl/input/smf.tsv \
      --genes VPREB3 \
      --chromosomes 22 \
      --celltypes B_intermediate \
      --test-type 'set'
```

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (1,000 individuals) - already generated
* Phase 1: OneK1K + BioHEART (2,000 individuals) - WGS done, scRNA-seq in progress
* Phase 2: TBC
