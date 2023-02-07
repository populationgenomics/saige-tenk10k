import logging

import hail as hl
from cpg_utils.hail_batch import (
    # copy_common_env,
    # dataset_path,
    # get_config,
    init_batch,
    output_path,
    # remote_tmpdir,
)

# similar to https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py
def filter_variants(
    mt_path: str,  # "mt/v7.mt"
    samples: list[str],
    output_mt_path: str,  # "tob_wgs_rv/densified_rv_only.mt"
    grm_plink_file: str,
):
    """Subset hail matrix table

    Input:
    joint call hail matrix table
    set of samples for which we have scRNA-seq data

    Output:
    subset hail matrix table, containing only variants that:
    1) are not ref-only, 2) biallelic, 3) meet QC filters, 4) are rare (MAF<5%)

    also, plink file containing variants that satisfy 1),2),3)
    but that are common (MAF>1%) to build sparse GRM
    """
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # subset to relevant samples (samples we have scRNA-seq data for)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # densify
    mt = hl.experimental.densify(mt)

    # filter out low quality variants and consider biallelic SNPss only (no multi-allelic, no ref-only, no indels)
    mt = mt.filter_rows(  # check these filters!
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNPs
    )

    mt = hl.variant_qc(mt)
    # filter common (enough) variants to build sparse GRM
    grm_mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > 0.01) & (mt.variant_qc.AF[1] < 1)
        | (mt.variant_qc.AF[1] < 0.99) & (mt.variant_qc.AF[1] > 0)
    )
    logging.info(f"Number of variants exported to build GRM: {grm_mt.count()[0]}")

    # export to plink common variants only for sparse GRM
    from hail.methods import export_plink

    export_plink(grm_mt, grm_plink_file, ind_id=grm_mt.s)

    # filter rare variants only (MAF < 5%)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.AF[1] > 0)
        | (mt.variant_qc.AF[1] > 0.95) & (mt.variant_qc.AF[1] < 1)
    )
    mt.write(output_mt_path, overwrite=True)
    logging.info(f"Number of rare (freq<5%) and QCed biallelic SNPs: {mt.count()[0]}")


# same as https://github.com/populationgenomics/cellregmap-pipeline/blob/main/batch.py
def get_promoter_variants(
    mt_path: str,  # output path from function above
    ht_path: str,  # open chromatin annos in same file??
    gene_details: dict[str, str],  # output of make_gene_loc_dict
    window_size: int,
    plink_file: str,  # "tob_wgs_rv/pseudobulk_rv_association/plink_files/GENE"
):
    """Subset hail matrix table

    Input:
    mt_path: path to already subsetted hail matrix table
    ht_path: path to VEP HT
    gene_details: dict of info for current gene
    window_size: int, size of flanking region around genes
    plink_file: str, file prefix for writing plink data

    Output:
    For retained variants, that are: 1) in promoter regions and
    2) within 50kb up or down-stream of the gene body (or in the gene body itself)
    (on top of all filters done above)

    returns nothing
    """

    # read hail matrix table object (pre-filtered)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    gene_name = gene_details["gene_name"]

    # get relevant chromosome
    chrom = gene_details["chr"]

    # subset to window
    # get gene body position (start and end) and build interval
    left_boundary = max(1, int(gene_details["start"]) - window_size)
    right_boundary = min(
        int(gene_details["end"]) + window_size,
        hl.get_reference("GRCh38").lengths[chrom],
    )
    # get gene-specific genomic interval
    gene_interval = f"{chrom}:{left_boundary}-{right_boundary}"
    logging.info(f"Interval considered: {gene_interval}")  # "chr22:23219960-23348287"

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome="GRCh38")]
    )
    mt_path = output_path(f"{gene_name}_in_window.mt", "tmp")
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f"Number of variants within interval: {mt.count()[0]}")

    # annotate using VEP
    vep_ht = hl.read_table(ht_path)
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)

    # filter variants found to be in promoter regions
    mt = mt.filter_rows(
        mt.vep.regulatory_feature_consequences["biotype"].contains("promoter")
    )
    promoter_path = output_path(f"{gene_name}_promoter_variants.mt", "tmp")
    mt = mt.checkpoint(promoter_path, overwrite=True)  # checkpoint
    logging.info(
        f"Number of rare (freq<5%) QC-passing, biallelic SNPs in promoter regions: {mt.count()[0]}"
    )

    # add aditional filtering for all regions that are open
    # these flags will need to be different for each cell type
    # but can be dealth with upon running (setting annotations in a way that only open regions for that chromatin are run)

    # export this as a Hail table for downstream analysis
    # consider exporting the whole MT?
    ht_path = output_path(
        f"summary_hts/{gene_name}_rare_promoter_summary.ht", "analysis"
    )
    ht = mt.rows()
    ht.write(ht_path, overwrite=True)

    # export MT object to PLINK (promoter variants)
    # pylint: disable=import-outside-toplevel
    from hail.methods import export_plink

    export_plink(mt, plink_file, ind_id=mt.s)
