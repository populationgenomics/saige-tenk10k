# Hail batch workflow to run SAIGE-QTL on TenK10K data

This is a Hail batch pipeline to run [SAIGE-QTL](https://github.com/weizhou0/qtl) on CPG's GCP, to map associations between both common and rare genetic variants and single-cell gene expression from peripheral mononuclear blood cells (PBMCs).
First, this will be run on the TOB (AKA OneK1K) and then BioHEART datasets as part of phase 1 of the TenK10K project (see *Data* below), but in the future all datasets within OurDNA (with scRNA-seq + WGS data) will go through this pipeline as well.

The pipeline is split into three parts, to make for more flexible usage:

1. Genotype processing (SNV): this involves variant QC and selection of the WGS data, and genotype file preparation specifically for common and rare single-nucleotide variants and indels (VCF files), as well as plink files for only a subset of 2,000 variants that is used for some approximations within SAIGE-QTL
2. Expression (phenotype) processing: this involves processing of the scRNA-seq data, inclusion of covariates, and preparation of the pheno_cov files (one per gene, cell type) and cis window files (one per gene)
3. Association testing: prepare and run SAIGE-QTL commands for association mapping using inputs generated in the first two parts.

Additionally, two helper scripts are also part of this pipeline:

* one to extract sample covariates which feeds into the expression processing script where they are combined with expression-based covariates
* one to make group files that are necessary for rare variant association testing

## Genotypes preprocessing

Script: get_genotype_vcf.py

Variant selection for VCF files:

* variants that are: i) QC-passing, ii) not ref-ref variants, and iii) not indels or multi-allelic SNPs (when run with --exclude-indels and --exclude-multiallelic).
* variants that are common at a set threshold (MAF > T) in our population (by default, T=0.01)
* one per chromosome

Variant selection for PLINK files for variance ratio estimation (VRE):

* variants that are: i) QC-passing, ii) not ref-ref variants, and iii) not indels or multi-allelic SNPs.
* variants that are not rare (MAC > 20) in our population
* random subset of 2,000 variants across all chromosomes

Inputs:

* joint call VDS object (TOB + BioHEART) after variant and sample QC has been applied.
* (for now, also gets tables of related individuals to exclude)

Outputs:

* VCF files containing all retained common variants (one per chromosome) + corresponding index files (`.csi`)
* VCF files containing all retained rare variants (one per chromosome) + corresponding index files (`.csi`)
* plink object for only 2,000 variants (minor allele count>20), after LD pruning - this is for the estimation of the variance ratio (VR plinks)

Notes: SAIGE-QTL allows numeric chromosomes only, so both the bim and the VCF files are modified in this script to remove the 'chr' notation (so that e.g. 'chr1' becomes '1').

## Get sample covariates

Script: `get_sample_covariates.py`

## Gene expression preprocessing

Script: `get_anndata.py`

Inputs:

* scanpy AnnData object (one per chromosome and cell type, TOB + BioHEART)
* cell covariate file (one per cell type, TOB + BioHEART)
* sample covariate file generated in get_samples_covariates.py

Outputs:

* TSV phenotype covariate files (one per gene, cell type)
* TSV gene cis window file (one per gene)

Notes: as before, we remove 'chr' from the chromosome name in the gene cis window file.
Additionally, we turn hyphens ('-') into underscores ('_') in the gene names.
Both the AnnData objects and cell covariate files are generated on Garvan's HPC and copied over to GCP.

## Make group file

Script: `make_group_file.py`

Inputs:

* rare variant VCF file
* cis window file

Outputs

* group files (one per gene)


## SAIGE-QTL association pipeline

Script: `saige_assoc.py`

Inputs:

* PLINK genotype files for VRE estimation (one only)
* VCF genotype files for SNP testing (one per chromosome, common variants)
* TSV phenotype covariate files for expression + covariate info (one per gene + cell type combination)
* TSV gene cis window file to specify what genomic region to test (one per gene)

Outputs:

* single-variant raw p-values (one per gene, cell type)
* association summary statistics (gene-corrected p-values summarised, one per cell type)

## SAIGE-QTL RV association pipeline

Script: `saige_assoc_set_test.py`

Inputs:

* PLINK genotype files for VRE estimation (one only)
* VCF genotype files for SNP testing (one per chromosome, rare variants)
* TSV phenotype covariate files for expression + covariate info (one per gene + cell type combination)
* TSV marker group file to specify what variants to test, and what weights and annotations to use (one per gene)

Outputs:

* set-based raw p-values (one per gene, cell type)
* association summary statistics (gene-corrected p-values summarised, one per cell type)

### SAIGE-QTL parameters explanation

Clarifying the reasoning behind the parameters / flags used to run SAIGE-QTL.
Most of these are (or will be) included in the official [documentation](https://weizhou0.github.io/SAIGE-QTL-doc/).

Note: some of these are provided as arguments in the scripts (`saige_assoc.py`, `saige_assoc_set_test.py`), but most are provided as a separate config file (`saige_assoc_test.toml`).

Fit null model ([step 1](https://weizhou0.github.io/SAIGE-QTL-doc/docs/step1.html)):

* ```pheno_file```: path specifying the location of the phenotype covariate file described above (build during part 2 of the pipeline)
* ```cov_col_list```: string specifying the columns of the pheno_file that should be used as covariates (separated by a comma, no spaces)
* ```sample_cov_col_list```: same as above, but specifying only, out of the columns above, the ones that are well defined at individual level (e.g., sex, age, ancestry PCs). Both this and the above need to be specified, and this is always a subset of the above, which allows individual-level covariates to be processed more cheaply.
* ```sample_id_pheno```: specify the column that represents individual IDs
* ```output_prefix```: path to where the output files from step 1 (which will be read by step 2) should be written to
* ```plink_path```: path to VRE plink files (specify just the prefix, but a .bim, .fam, and .bed files with the same prefix and in the same location should exist -- these are built in part 1)
* ```pheno_col```: specify the column that should be used as phenotype, in our case the gene to test
* ```trait_type```: specify the model to be used, ```count``` is the Poisson model which should be used here.
* ```skip_vre```: boolean specifying whether the variance ratio estimation should be run or skipped, should always be false (Note that because the syntax is different between Python and R we encode this as the string ```FALSE``` instead of the boolean ```False```)
* ```pheno_remove_zeros```: option to remove 0s from phenotype vector (default: ```FALSE``` as it does not make sense for the very sparse scRNA-seq data)
* ```use_sparse_grm_null```: use sparse GRM to account for relatedness. This is implemented but would require a step0 in the pipeline to first construct this, which is not there at the moment (default: ```FALSE```)
* ```use_grm_null```: same as above, but without the sparse option (default: ```FALSE```)
* ```is_cov_offset```: if there are no covs, add an offset of ones (never really the case for us, default: ```FALSE```)
* ```is_cov_transform```: transform (explain) covariates (default: ```TRUE```)
* ```skip_model_fitting```: boolean (default: ```FALSE```)
* ```tol```: convergence tolerance (default 0.00001, which works well in our hands)
* ```maxiterPCG```: convergence max number of iterations (default 500 but increase to 5,000 if job does not converge)
* ```is_overwrite_vre_file```: if the file already exists, skip or overwrite, default is the latter (default: ```TRUE```)

Single-variant association testing ([common variants step 2](https://weizhou0.github.io/SAIGE-QTL-doc/docs/single_step2.html)):

* ```vcf_file```: path to VCF file containing genetic variants to be tested
* ```vcf_file_index```: corresponding .csi index file (not .tbi)
* ```vcf_field```: DS for dosages, GT for genotypes (default = 'GT')
* ```saige_output_file```: path to output file
* ```chrom```: chromosome to be tested
* ```cis_window_file```: path to file specifying cis window / region to test (generated in part 2 of the pipeline, get anndata script)
* ```gmmat_model_path```: path to estimated null model (.rda) generated in step 2
* ```variance_ratio_path```: path to variance ratio txt file generated in step 1
* ```min_maf```: minimum minor allele frequency (MAF) (default: 0)
* ```min_mac```: minimum minor allele count (MAC) (default: 5)
* ```loco_bool```: boolean specifying whether leave-one-chromosome-out should be used (default: ```FALSE```)
* ```n_markers```: internal parameter to batchify variants tested (default: 10000),
* ```spa_cutoff```: internal parameter to do with the saddlepoint approximation, does not make much of a difference for us (default: 10000).

Obtain gene-level p-values ([common variants only, step 3](https://weizhou0.github.io/SAIGE-QTL-doc/docs/gene_step3.html))

* ```gene_name```: gene to aggregate values for,
* ```saige_sv_output_file```: path to output from step 2 (input here),
* ```saige_gene_pval_output_file```: path to output (step 3).

Set-based association testing ([rare variants step 2](https://weizhou0.github.io/SAIGE-QTL-doc/docs/set_step2.html)):

* ```vcf_file```: path to VCF file containing genetic variants to be tested
* ```vcf_file_index```: corresponding .csi index file (not .tbi)
* ```vcf_field```: DS for dosages, GT for genotypes (default = 'GT')
* ```saige_output_file```: path to output file (if kept the same as single-variant test, step 1 only needs to be run once)
* ```chrom```: chromosome to be tested
* ```group_file```: path to file specifying variants to test (generated in the make_group_file.py script)
* ```gmmat_model_path```: path to estimated null model (.rda) generated in step 2
* ```variance_ratio_path```: path to variance ratio txt file generated in step 1
* ```min_maf```: minimum minor allele frequency (MAF) (default: 0)
* ```min_mac```: minimum minor allele count (MAC) (default: 5)
* ```loco_bool```: boolean specifying whether leave-one-chromosome-out should be used (default: ```FALSE```)
* ```n_markers```: internal parameter to batchify variants tested (default: 10000),
* ```spa_cutoff```: internal parameter to do with the saddlepoint approximation, does not make much of a difference for us (default: 10000),

## To run

Instructions to run each component of the pipeline using analysis runner are provided at the top of each script.

## Data

TenK10K is matched single-cell RNA-seq (scRNA-seq) and whole-genome sequencing (WGS) data from up to 10,000 individuals:

* Phase 0: OneK1K data only (~1,000 individuals) - already generated (both WGS and scRNA-seq, though with an older technology)
* Phase 1: OneK1K + BioHEART (~2,000 individuals) - WGS done, scRNA-seq in progress (almost done)
* Phase 2/final: aim is ~ 10,000 individuals from the (extended) TOB/OneK1K, BioHEART, ADAPT, LBIO and AIM cohorts (nothing generated besides current stage of Phase 1)

## Additional resources

* [SAIGE-QTL preprint](https://www.medrxiv.org/content/10.1101/2024.05.15.24307317v1)
* [OneK1K paper](https://www.science.org/doi/10.1126/science.abf3041)
* [SAIGE-QTL pipeline flowchart GSlides](https://docs.google.com/presentation/d/1OhNiA6DaP9ZGlAbh8uZuZWzvrrr_QwvJwJ_lVPBZoic/edit#slide=id.g25daf727307_0_102)
* [SAIGE-QTL pipeline notes GDoc](https://docs.google.com/document/d/1t11VafeU1THA4X58keHd5oPVglTYiY3DKC7P05GHCzw/edit)
