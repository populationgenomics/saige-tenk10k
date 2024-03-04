# Hail batch workflow to run SAIGE-QTL on TenK10K data

This is a hail batch pipeline to run the new [QTL version of SAIGE](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between both common and rare genetic variants and single-cell gene expression from peripheral mononuclear blood cells (PBMCs).
First, this will be run on the TOB (AKA OneK1K) and then BioHEART datasets as part of phase 1 of the TenK10K project (see *Data* below), but in the future all datasets within OurDNA (with scRNA-seq + WGS data) will go through this pipeline as well.

The pipeline is split into three parts, to make for more flexible usage:

1. Genotype processing (SNV): this involves variant QC and selection of the WGS data, and genotype file preparation specifically for common and rare single-nucleotide variants (VCF files), as well as plink files for only a subset of 2,000 variants that is used for some approximations within SAIGE-QTL
2. Expression (phenotype) processing: this involves processing of the scRNA-seq data, inclusion of covariates, and preparation of the pheno_cov files (one per gene, cell type) and cis window files (one per gene)
3. Association testing: prepare and run SAIGE-QTL commands for association mapping using inputs generated in the first two parts.

## Genotypes preprocessing

Script: get_genotype_vcf.py

Variant selection for VCF files:

* variants that are: i) QC-passing, ii) not ref-ref variants, and iii) not indels or multi-allelic SNPs (when run with --exclude-indels and --exclude-multiallelic).
* variants that are common (MAF > 0.01) in our population
* one per chromosome

Variant selection for PLINK files for variance ratio estimation (VRE):

* variants that are: i) QC-passing, ii) not ref-ref variants, and iii) not indels or multi-allelic SNPs.
* variants that are not rare (MAC > 20) in our population
* random subset of 2,000 variants across all chromosomes

Inputs:

* joint call VDS object (TOB + BioHEART)

Outputs:

* VCF file containing all retained common variants (one per chromosome) + corresponding index file (.csi)
* plink object for only 2,000 variants (minor allele count>20), after LD pruning - this is for the estimation of the variance ratio (VR plinks)

## Gene expresion preprocessing

Script: get_anndata.py

Inputs:

* scanpy AnnData object (one per chromosome and cell type, TOB + BioHEART)

Outputs:

* TSV phenotype covariate files (one per gene, cell type)
* TSV gene cis window file (one per gene)

## SAIGE-QTL association pipeline

Script: saige_assoc.py

Inputs:

* PLINK genotype files for VRE estimation (one only)
* VCF genotype files for SNP testing (one per chromosome)
* TSV phenotype covariate files for expression + covariate info (one per gene + cell type combination)
* TSV gene cis window file to specify what variants to test (one per gene)

Outputs:

* association summary statistics

### SAIGE-QTL parameters explanation

Fit null model ([step 1](https://weizhou0.github.io/SAIGE-QTL-doc/docs/step1.html)):

* ```pheno_file```: path specifying the location of the phenotype covariate file described above (build during part 2 of the pipeline) 
* ```cov_col_list```: string specifying the columns of the pheno_file that should be used as covariates (separated by a comma)
* ```sample_cov_col_list```: same as above, but specifying only, out of the columns above, the ones that are well defined at individual level (e.g., sex, age)
* ```sample_id_pheno```: specify the column that represents individual IDs
* ```output_prefix```: path to where the output files from step 1 (which will be read by step 2) should be written to
* ```plink_path```: path to VRE plink files (specify just the prefix, but a .bim, .fam, and .bed files with the same prefix and in the same location should exist (these are built in part 1)
* pheno_col: str = 'y',
* trait_type: str = 'count',
* skip_vre: str = 'FALSE',
* pheno_remove_zeros: str = 'FALSE',
* use_sparse_grm_null: str = 'FALSE',
* use_grm_null: str = 'FALSE',
* is_cov_offset: str = 'FALSE',
* is_cov_transform: str = 'TRUE',
* skip_model_fitting: str = 'FALSE',
* tol: float = 0.00001,
* is_overwrite_vre_file: str = 'TRUE',


## To run

Instructions to run each component of the pipeline using analysis runner are provided at the top of each script

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (~1,000 individuals) - already generated (both WGS and scRNA-seq, though with an older technology)
* Phase 1: OneK1K + BioHEART (~2,000 individuals) - WGS done, scRNA-seq in progress (almost done)
* Phase 2/final: aim is ~ 10,000 individuals from the (extended) TOB/OneK1K, BioHEART, ADAPT, LBIO and AIM cohorts (nothing generated besides current stage of Phase 1)

## Additional resources

* [SAIGE-QTL pipeline flowchart GSlides](https://docs.google.com/presentation/d/1OhNiA6DaP9ZGlAbh8uZuZWzvrrr_QwvJwJ_lVPBZoic/edit#slide=id.g25daf727307_0_102)
* [SAIGE-QTL pipeline notes GDoc](https://docs.google.com/document/d/1t11VafeU1THA4X58keHd5oPVglTYiY3DKC7P05GHCzw/edit)
