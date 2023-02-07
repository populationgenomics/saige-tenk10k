
from cpg_utils.hail_batch import (
    get_config,
    remote_tmpdir,
)

import hailtop.batch as hb

HAIL_IMAGE = 'blablabla'

sb = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )

batch = hb.Batch('SAIGE QTL pipeline', backend=sb)

# FIRST THREE STEPS DONE ONLY ONCE

# input: all variants called by WGS (MT)
# output: filtered MT i) QC, ii) no ref-ref
filter_job = batch.new_python_job(name='MT filter job')
filter_job.image(HAIL_IMAGE)  # what image??
filter_job.call(params)

# input: filtered MT
# output: filtered MT common variants only
grm_variants_job = batch.new_python_job(name='GRM variants filter job')
grm_variants_job.depends_on(filter_job)
grm_variants_job.image(HAIL_IMAGE)
grm_variants_job.call(params)

# input: filtered MT
# output: filtered MT rare variants only
rv_filter_job = batch.new_python_job(name='rare variants filter job')
rv_filter_job.depends_on(filter_job)
rv_filter_job.image(HAIL_IMAGE)
rv_filter_job.call(params)
