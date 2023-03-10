# Hail batch workflow to run SAIGE-QTL on TenK10K data

This is a hail batch pipeline to run the new [QTL version of SAIGE](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between common and rare genetic variants and single-cell gene expression from blood.
First, this will be run on the TOB and then BioHEART datasets as part of the TenK10K project (see *Data* below), but in the future all datasets within OurDNA (with scRNA-seq + WGS data) will go through this pipeline as well.

* **Plan A**: at present, just adapting our [CellRegMap Hail batch pipeline](https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py) hoping it can run R / external code smoothly.
* **Plan A**: at present, just adapting our [CellRegMap Hail batch pipeline](https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py) hoping it can run R / external code smoothly.
* **Plan B**: if that fails, we may need to adapt [Konrad K's UKBB exomes analysis github](https://github.com/Nealelab/ukb_exomes), underlying [this paper](https://www.sciencedirect.com/science/article/pii/S2666979X22001100), or at least using [these python wrappers for SAIGE](https://github.com/Nealelab/ukb_common/blob/master/utils/saige_pipeline.py).

## Plan A

### Genotypes preprocessing (once per cohort)

Function name: filter_variants.

Hail query to filter WGS object to i) QC-passing, ii) non ref-ref variants, and considering only samples with scRNA-seq data.
At present, this also filters out indels and multi-allelic SNPs, which we may want to relax.

It outputs three objects:

* MT object, rare (freq<5%) variants
* MT object, common (freq>1%) variants
* plink object for only 2,000 variants (MAC>20), after LD pruning - this is for the estimation of the variance ratio (VR plinks)

<!-- # skip for now - unrelated individuals
* SAIGE R script to create sparse GRM
  * just once for all individuals, all variants after LD-pruning, and MAF>1% -->

### Expression preprocessing (once per cell type, chromosome)

Function name: prepare_pheno_cov_file.

Python script to combine expression (pheno), covariates into a single pheno_cov file as input to saige-qtl.

Inputs:

* chromosome and cell type-specific scanpy (AnnData) objects, single-cell expression for all cells, all genes for that chromosome and all cells assigned to that cell type (sctransformed sc-counts)
* covariates: combination of 1) exploded donor-level covariates (e.g. sex) and 2) cell-level covariates (PCs, batch)
* sample mapping file: matching donor ID (onek1k vs CPG) and cell barcodes, including cell type labels

Output file (one per cell type):

* text file concatenating covs, expression of genes, and individual id (same as plink files) for all cells from that cell type

i.e.,

ind_id | cov1 | ... | covM | gene1 | ... | geneN

Note: if these files get too large, create one pheno_cov file per gene instead, i.e.,

File1:

ind_id | cov1 | ... | covM | gene1

FileN:

ind_id | cov1 | ... | covM | geneN

### Variant selection (once per gene)

Function name: get_variants_to_test.

Hail query to filter object to relevant variants (for a specific gene)

* for common variants, only genomic proximity (+/-100kb)
* for rare variants, i) genomic proximity (+/-50kb), 2) regulatory consequences (vep), 3) open chromatin (any cell type)

Outputs:

* plink files (.bed, .bim, .fam) for common variants
* plink files (.bed, .bim, .fam) for rare variants
* group file for rare variants. For each region(=gene?)
  * row1: region1 | var    | 1:100010:A:C | 1:100016:A:C
  * row2: region1 | anno   | promoter     | enhancer,promoter
  * row3: region1 | weight | 0.5          | 0.3

Note: weights are not necessary, annotations can be multiple, separated by comma.

### Run association (for each gene, cell type combination)

Function names: build_fit_null_command, build_run_association_test_command.

This sets up R commands to actually run SAIGE-QTL, which is run in two steps.

Step1 fits a null model (no variants yet) for the gene of interest, given the appropriate covariates, whereas step2 actually runs associations for the variants provided.
Step 2 comes in two flavours, using either a region-based test (best suited for testing rare variants' effects), or a single-variant tests (typically used for common variants).

The command runs an R script, but otherwise behaves like a command line tool.
As such, for both steps 1 and 2, we consider,

  * Python function to build command (which just adds to a string, i.e., "Rscript myscript.R --arg1 par1 --arg2 par2..")
  * Job to actually run the command with the appropriate arguments

#### Step1: Fit null model

Function name: build_fit_null_command.

Inputs:

* Pheno cov file
* VR plinks
* info on what gene, what covs, what likelihood to use

Outputs:

* model object (.rda)
* variance ratio estimate (.txt)
* optional: log

#### Step2: Run association

Function name: build_run_set_test_command.

Inputs:

* model object (.rda)
* variance ratio estimate (.txt)
* genotypes to test (plink files: .bed, .bim, .fam)
* 'set-test' only: group file with group definitions and annotations

Output:

* txt file with all results
* optional: log

### Results aggregation (once per cell type)

Function name: summarise_association_results

Python to aggregate all results from step 2 above.
Result tables are gene-specific, this aggregated them all both to have a single-summary and to appropriately perform multiple testing correction.
If both a set-test and a single-variant association test have been performed, this will results in 2 X n_cell_types final result tables.

Input:

* cell type of interest
* (cell type and) gene-level specific result tables (.txt)

Output:

* cell type-specific summary statistics across all genes (.txt)

### To run (this may need to be updated)

```bash
analysis-runner \
    --dataset tob-wgs \
    --access-level standard \
    --output-dir 'tob_wgs_genetics/saige_qtl/' \
    --image 'australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:dev' \
    --description 'Saige QTL batch job' \
    python3 saige_tenk10k.py \
      --input-files-prefix 'tob_wgs_genetics/saige_qtl/input' \
      --sample-mapping-file-tsv 'tob_wgs_genetics/saige_qtl/input/smf.tsv' \
      --genes 'VPREB3' \
      --chromosomes 22 \
      --celltypes 'B_intermediate' \
      --test-type 'set'
```

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (1,000 individuals) - already generated
* Phase 1: OneK1K + BioHEART (2,000 individuals) - WGS done, scRNA-seq in progress (new kit, which should result in many more cells per individual)
* Phase 2/final: aim is ~ 10,000 individuals from the (extended) TOB/OneK1K, BioHEART, ADAPT, LBIO and AIM cohorts (nothing generated besides current stage of Phase1)
