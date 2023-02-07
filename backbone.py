
from cpg_utils.hail_batch import (
    get_config,
    remote_tmpdir,
)

import hailtop.batch as hb

HAIL_IMAGE = 'blablabla'  # needs hail query
SAIGE_QTL_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/saige-qtl'
PY_IMAGE = 'blabla'  # does not, but may need scanpy (or Seurat??)

sb = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )

batch = hb.Batch('SAIGE QTL pipeline', backend=sb)

# FIRST STEPS DONE ONLY ONCE

# input: all variants called by WGS (MT)
# output: filtered MT i) QC, ii) no ref-ref, iii) scRNA-seq individuals only
filter_job = batch.new_python_job(name='MT filter job')
filter_job.image(HAIL_IMAGE)  # what image??
filter_job.call(params)

# input: filtered MT
# output: GRM MT i.e., filtered MT common variants only?
grm_variants_job = batch.new_python_job(name='GRM variants filter job')
grm_variants_job.depends_on(filter_job)
grm_variants_job.image(HAIL_IMAGE)
grm_variants_job.call(params)

# input: GRM MT
# output: Sparse GRM object
sparse_grm_job = batch.new_job(name='Create sparse GRM')
sparse_grm_job.depends_on(grm_variants_job)
sparse_grm_job.image(SAIGE_QTL_IMAGE)
sparse_grm_job.call(params)

# input: filtered MT
# output: filtered MT rare variants only
rv_filter_job = batch.new_python_job(name='rare variants filter job')
rv_filter_job.depends_on(filter_job)
rv_filter_job.image(HAIL_IMAGE)
rv_filter_job.call(params)

# input: GRM MT or filtered MT or rare variants(check with Wei)
# output: plink object for a subset of variants for variance ratio estimation
vr_geno_job = batch.new_python_job(name='Variance ratio subset job')
vr_geno_job.depends_on(filter_job)
# vr_geno_job.depends_on(grm_variants_job)
# vr_geno_job.depends_on(rv_filter_job)
vr_geno_job.image(HAIL_IMAGE)
vr_geno_job.call(params)


# identify relevant genes,
# i.e., expressed in (=present in pheno file for)
# at least one cell type
# or provided as inputs
relevant_genes = ['gene1', 'gene2']

# ONCE PER GENE, ACROSS ALL CELL TYPES

# input: RV filtered MT, gene, optional: regulatory category
# output: MT i) filtered, ii) rare + iii) proximal, iv) regulatory (possibly only some category)
for gene in relevant_genes:
    gene_job = batch.new_python_job(f'Extract proximal regulatory variants for: {gene}')
    gene_job.depends_on(filter_job)
    gene_job.image(HAIL_IMAGE)
    gene_job.call(params)


# ONCE FOR EVERY GENE AND CELL TYPE COMBINATION
celltypes = ['celltype1', 'celltype2']

# jobs below may be combined but add them separately at this stage for clarity
for celltype in celltypes:
    # obtain cell type expressed genes
    for gene in celltype_genes:

        # input: gene ID, phenotype file, covariate file
        # output: pheno_cov file
        pheno_cov_job = batch.new_python_job(f'Make pheno cov file for: {gene}, {celltype}')
        # does not depend on any other jobs
        pheno_cov_job.image(PY_IMAGE)
        pheno_cov_job.call(params)

        # input: sparse GRM, pheno_cov file, subset plink files
        # output: null model object, variance ratio (VR) estimate file
        fit_null_job = batch.new_python_job()
        # syntax below probably does not work
        dependencies = [sparse_grm_job, vr_geno_job, pheno_cov_job]
        fit_null_job.depends_on(dependencies)
        fit_null_job.image(SAIGE_QTL_IMAGE)
        # python job creating Rscript command
        fit_null_job.call(
            get_step1_command,
            params)
        # regular job submitting the Rscript command to bash
        fit_null_job.call(
            run_saige_step1,
            params)

        # input: null model object, VR file, gene specific genotypes
        run_association_job = batch.new_python_job()
        dependencies = [fit_null_job, gene_job]
        run_association_job.depends_on(dependencies)
        run_association_job.image(SAIGE_QTL_IMAGE)
        # python job creating Rscript command
        run_association_job.call(
            get_step2_command,
            params)
        # regular job submitting the Rscript command to bash
        run_association_job.call(
            run_saige_step2,
            params)
