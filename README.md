# Hail batch workflow to run SAIGE-QTL on TenK10K data

This is a Hail [batch](https://hail.is/docs/batch/index.html) pipeline to run [SAIGE-QTL](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between both common and rare genetic variants and single-cell gene expression from peripheral mononuclear blood cells (PBMCs).
First, this will be run on the TOB (AKA OneK1K) and then BioHEART datasets as part of phase 1 of the TenK10K project (see *Data* below), but in the future all datasets within OurDNA (with scRNA-seq + WGS data) will go through this pipeline as well.

The pipeline is split into three main parts, to make for more flexible usage:

1. Genotype processing (SNVs and indels): this involves sample and variant QC of the WGS data, and genotype file preparation specifically for common and rare single-nucleotide variants and indels (VCF files), as well as plink files for only a subset of 2,000 variants that is used for some approximations within SAIGE-QTL
2. Expression (phenotype) processing: this involves processing of the scRNA-seq data, inclusion of covariates, and preparation of the phenotype + covariate files (one per gene, cell type) and cis window files (one per gene)
3. Association testing: prepare and run SAIGE-QTL commands for association mapping using inputs generated in the first two parts.

Additionally, two helper scripts are also part of this pipeline:

* one to extract sample covariates which feed into the expression processing script where they are combined with expression-based covariates
* one to make group files that are necessary for rare variant association testing

## Genotypes preprocessing

Script: [`get_genotype_vcf.py`](get_genotype_vcf.py)

Variant selection for VCF files:

* variants that are: i) QC-passing, ii) not ref-ref variants, and iii) not indels or multi-allelic SNPs (when run with `--exclude-indels` and `--exclude-multiallelic`).
* separately, variants that are common or rare at a set threshold (MAF >= T, MAF < T, respectively) in our population (by default, T=0.01)
* two per chromosome (one for common, one for rare variants)

Variant selection for PLINK files for variance ratio estimation (VRE):

* variants that are: i) QC-passing, ii) not ref-ref variants, and iii) not indels or multi-allelic SNPs.
* variants that are not rare (MAC > 20) in our population
* random subset of 2,000 variants across all chromosomes

Inputs:

* joint call VDS object (TOB + BioHEART) after variant and sample QC has been applied (ideally, not yet).
* (for now, also gets tables of related individuals to exclude)

Outputs:

* VCF files containing all retained common variants (one per chromosome) + corresponding index files (`.csi`)
* VCF files containing all retained rare variants (one per chromosome) + corresponding index files (`.csi`)
* plink objects (`.bed`, `.bim`, `.fam`) for only 2,000 variants (minor allele count > 20), after LD pruning - this is for the estimation of the variance ratio (VRE plink files)

### Notes (get_genotype_vcf.py)

SAIGE-QTL allows numeric chromosomes only, so both the `.bim` and the VCF files are modified in this script to remove the 'chr' notation (so that e.g., 'chr1' becomes '1').

For the VRE estimation, we need to select 2,000 (this number can be customised) variables at random, except we need them to be common enough (MAC>20, also customisable), and we want them to be somewhat independent, to be more representative of the entire cohort.
We also subset to only variants on autosome (1-22) chromosomes.
To achieve this, we perform LD pruning of the variables so that variants in strong LD get pruned to one in each LD block.
Because the LD pruning operation is very costly, we downsample the variants first (to 1% by default).
The LD pruning parameters can be user-adjusted, with default values as described in the [methods's documentation](https://hail.is/docs/0.2/guides/genetics.html#pruning-variants-in-linkage-disequilibrium).

## Get sample covariates

Script: [`get_sample_covariates.py`](get_sample_covariates.py)

Inputs:

* sex info for the cohort(s) of interest
* genotype principal components for the cohort(s) of interest
* (information from metamist, at the moment just age)

Outputs:

* TSV sample covariate file (one per cohort)

### Notes (get_sample_covariates.py)

There is the option to fill in missing values for sex (0 for unknown, where 1 is male, 2 is female) and age (average age across the cohort).
Additionally, add a user-specified (default: 10) number of permuted IDs, where the individual ID is permuted at random, to assess calibration (by shuffling the individual IDs we disrupt any real association between genotype and phenotype, so we expect no significant associations left when testing).

## Gene expression preprocessing

Script: [`get_anndata.py`](get_anndata.py)

Inputs:

* scanpy AnnData object (one per chromosome and cell type, TOB + BioHEART)
* cell covariate file (one per cell type, TOB + BioHEART)
* sample covariate file generated in get_samples_covariates.py

Outputs:

* TSV phenotype covariate files (one per gene, cell type)
* TSV gene cis window file (one per gene)

### Notes (get_anndata.py)

As in `get_genotype_vcf.py`, we remove 'chr' from the chromosome name in the gene cis window file.
Additionally, we turn hyphens ('-') into underscores ('_') in the gene names.
Both the AnnData objects and cell covariate files are generated on Garvan's HPC and copied over to GCP.
A note that the `filter_lowly_expressed_genes` method will remove lowly-expressed genes that will not even get tested, which should be kept in mind when interpreting the results.

## Make group file

Script: [`make_group_file.py`](make_group_file.py)

Inputs:

* rare variant VCF file
* cis window file

Outputs

* group files (one per gene)

### Notes (make_group_file.py)

Option to include no weights or to compute weights that reflect the distance of each variant from the gene's transcription start site (`dTSS`).
Using one of the flags below it is possible to additionally test using equal weights.
We use no annotations for now (set to `null`).

## SAIGE-QTL common variant association pipeline

Script: [`saige_assoc.py`](saige_assoc.py)

Run this for single-variant tests (typically for common variants).

Inputs:

* PLINK genotype files for VRE estimation (one only)
* VCF genotype file (+ index file) for SNP testing (one per chromosome, common variants)
* TSV phenotype covariate files for expression + covariate info (one per gene + cell type combination)
* TSV gene cis window file to specify what genomic region to test (one per gene)

Outputs:

* single-variant raw p-values (one per gene, cell type)
* association summary statistics (ACAT gene-corrected p-values summarised, one per cell type)

## SAIGE-QTL rare variant association pipeline

Script: [`saige_assoc_set_test.py`](saige_assoc_set_test.py)

Run this for set-based tests (typically for rare variants).

Inputs:

* PLINK genotype files for VRE estimation (one only)
* VCF genotype files for SNP testing (one per chromosome, rare variants)
* TSV phenotype covariate files for expression + covariate info (one per gene + cell type combination)
* TSV marker group file to specify what variants to test, and what weights and annotations to use (one per gene)

Outputs:

* set-based test raw p-values (one per gene, cell type)
* if set to true, single-variant test raw p-values for all variants in the group also (one per gene, cell type)
* set-based association summary statistics (gene-level p-values summarised, one per cell type)

## SAIGE-QTL parameters explanation

Clarifying the reasoning behind the parameters / flags used to run SAIGE-QTL.
Most of these are (or will be) included in the official [documentation](https://weizhou0.github.io/SAIGE-QTL-doc/).

Note: some of these are provided as arguments in the scripts (`saige_assoc.py`, `saige_assoc_set_test.py`), but most are provided as a separate config file (`saige_assoc_test.toml`).

Fit null model ([step 1](https://weizhou0.github.io/SAIGE-QTL-doc/docs/step1.html)).

In script (`saige_assoc.py`, `saige_assoc_set_test.py`, using standard Snake Case naming convention as in the rest of the scripts):

* `pheno_file`: path specifying the location of the phenotype covariate file described above (built during part 2 of the pipeline).
* `output_prefix`: path to where the output files from step 1 (which will be read by step 2) should be written to.
* `plink_path`: path to VRE plink files (specify just the prefix, but a .bim, .fam, and .bed files with the same prefix and in the same location should exist -- these are built in part 1).
* `pheno_col`: specify the column that should be used as phenotype, in our case the gene to test.

In config (under `[saige.build_fit_null]` in `saige_assoc_test.toml`, using the Camel Case naming convention adopted in SAIGE-QTL):

* `covarColList`: string specifying the columns of the pheno_file that should be used as covariates (separated by a comma, no spaces).
* `sampleCovarColList`: same as above, but specifying only, out of the columns above, the ones that are well defined at individual level (e.g., sex, age, ancestry PCs). Both this and the above need to be specified, and this is always a subset of the above, which allows individual-level covariates to be processed more cheaply.
* `sampleIDColinphenoFile`: specify the column that represents individual IDs (note, for calibration analysis use one of the permuted ids here).
* `traitType`: specify the model to be used, `count` is the Poisson model which should be used here.
* `skipVarianceRatioEstimation`: boolean specifying whether the variance ratio estimation should be run or skipped, should always be false in our case (note that because the syntax of booleans is different between Python and R we encode this as the string `FALSE` instead of the boolean `False`).
* `isRemoveZerosinPheno`: option to remove 0s from phenotype vector (default: `FALSE` as it does not make sense for the very sparse scRNA-seq data, i.e. these are not missing data!).
* `useSparseGRMtoFitNULL`: use sparse GRM to account for relatedness. This is implemented but would require a step0 in the pipeline to first construct this, which is not there at the moment (default: `FALSE`).
* `useGRMtoFitNULL`: same as above, but without the sparse option (default: `FALSE`).
* `isCovariateOffset`: use covariates as offset allows to only estimate coefficients once and reduce cost (check that the model retains calibration, in which case set to `TRUE` for cost optimisation).
* `isCovariateTransform`: transform (explain) covariates (default: `TRUE`).
* `skipModelFitting`: boolean (default: `FALSE`).
* `tol`: convergence tolerance (default: 0.00001, which works well in our hands, we could test with a larger number to decrease cost).
* `maxiterPCG`: convergence max number of iterations (default: 500 but increase to 5,000 for the specific genes whose job does not converge).
* `IsOverwriteVarianceRatioFile`: if the file already exists, skip or overwrite, default is the latter (i.e. default: `TRUE`).

Single-variant association testing ([common variants step 2](https://weizhou0.github.io/SAIGE-QTL-doc/docs/single_step2.html)):

In script (`saige_assoc.py`):

* `vcf_file`: path to VCF file containing genetic variants to test. An index file that is called exactly the same, with `.csi` at the end, needs to exist at the same location, but is not passed in as an argument.
* `sv_output_path`: path to output file.
* `chrom`: chromosome to test.
* `cis_window_file`: path to file specifying cis window / region to test (generated in part 2 of the pipeline, get anndata script).
* `gmmat_model_path`: path to estimated null model (.rda) generated in step 1.
* `variance_ratio_path`: path to variance ratio txt file generated in step 1.

In config (under `[saige.sv_test]` in `saige_assoc_test.toml`):

* `vcfField`: DS for dosages, GT for genotypes (default: `GT`).
* `minMAF`: minimum minor allele frequency (MAF) (default: 0).
* `minMAC`: minimum minor allele count (MAC) (default: 5). Note that if this filter is discordant with the above one, the more stringent one will be applied (max between the two will be used).
* `LOCO`: boolean specifying whether leave-one-chromosome-out should be used (default: `FALSE`).
* `markers_per_chunk`: internal parameter to batchify variants tested (default: 10000). Confusingly the naming is Snake Case here.
* `SPAcutoff`: internal parameter to do with the saddlepoint approximation, does not make much of a difference for us (default: 10000).

Obtain gene-level p-values ([common variants only, step 3](https://weizhou0.github.io/SAIGE-QTL-doc/docs/gene_step3.html))

* `gene_name`: gene to aggregate values for,
* `saige_sv_output_file`: path to output from step 2 (input here),
* `saige_gene_pval_output_file`: path to output (step 3).

Set-based association testing ([rare variants step 2](https://weizhou0.github.io/SAIGE-QTL-doc/docs/set_step2.html)):

In script (`saige_assoc_set_test.py`):

* `vcf_file`: path to VCF file containing genetic variants to test. As above, an index file that is called exactly the same, with `.csi` at the end, needs to exist at the same location, but is not passed in as an argument.
* `set_output_path`: path to output file (if kept the same as single-variant test, step 1 only needs to be run once)
* `chrom`: chromosome to test.
* `group_file`: path to file specifying variants to test (generated in the `make_group_file.py` script).
* `gmmat_model_path`: path to estimated null model (`.rda`) generated in step 1.
* `variance_ratio_path`: path to variance ratio txt file generated in step 1.

In config (under `[saige.set_test]` in `saige_assoc_test.toml`):

* `vcfField`: DS for dosages, GT for genotypes (default: `GT`).
* `maxMAF_in_groupTest`: max MAF of variants that should be included in the set test (default: 0.1).
* `minMAF_in_groupTest_Exclude`: min MAF of variants to *exclude* from the set test (default: 0).
* `annotation_in_groupTest`: annotations from the group file to test as separate sets (default: `null`).
* `MACCutoff_to_CollapseUltraRare`: minor allele counts (MAC) cutoff under which variants are considered ultrarare and get collapsed (default: 10).
* `is_single_in_groupTest`: run single-variant tests for variants in group file as well (default: `TRUE`).
* `is_equal_weight_in_groupTest`: test using equal weights as well as using a Beta(1,25) distribution which is used by default; if set to true, both results get output (default: `TRUE`).

## To run

Instructions to run each component of the pipeline using analysis runner are provided at the top of each script.

Briefly, if one wanted to run both common and rare variant pipelines, the order of running would be:

1. get_genotype_vcf.py
2. get_sample_covariates.py (does not require any other part of the pipeline and can be run in parallel with 1)
3. get_anndata.py (requires 2)
4. saige_assoc.py (requires 1,3)
5. make_group_file.py (requires 1,3 and can be run in parallel with 4)
6. saige_assoc_set_test.py (requires 1,3,4 (so that step1 is only run once) and 5)

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (~1,000 individuals) - already generated (both WGS and scRNA-seq, though with an older technology)
* Phase 1: OneK1K + BioHEART (~2,000 individuals) - freeze 1 generated (~1,700 individuals), freeze 2 ongoing
* Phase 2/final: aim is ~ 10,000 individuals from the (extended) TOB/OneK1K, BioHEART, ADAPT, LBIO and AIM cohorts (nothing generated besides current stage of Phase 1)

## Additional resources

* [SAIGE-QTL preprint](https://www.medrxiv.org/content/10.1101/2024.05.15.24307317v1)
* [OneK1K paper](https://www.science.org/doi/10.1126/science.abf3041)
* [SAIGE-QTL pipeline flowchart GSlides](https://docs.google.com/presentation/d/1OhNiA6DaP9ZGlAbh8uZuZWzvrrr_QwvJwJ_lVPBZoic/edit#slide=id.g25daf727307_0_102)
* [SAIGE-QTL pipeline notes GDoc](https://docs.google.com/document/d/1t11VafeU1THA4X58keHd5oPVglTYiY3DKC7P05GHCzw/edit)

## Abbreviations

* ACAT: Aggregated Cauchy Association Test (PMID: 30849328)
* AIM: Australia IBD Microbiome study
* CPG: centre for population genomics, the centre at which this work is being done
* GCP: Google Cloud Platform
* HPC: high performance computing
* LBIO: liquid biopsies, a cancer cohort
* LD: linkage disequilibrium
* MAC: minor allele count
* MAF: minor allele frequency
* scRNA-seq: single-cell RNA sequencing
* SNP: single-nucleotide polymorphism
* TOB: Tasmanian Ophtalmic Biobank
* TSV: tab-separated file
* VCF: Variant Call Format
* VRE: variant ration estimate
* WGS: whole-genome sequencing
