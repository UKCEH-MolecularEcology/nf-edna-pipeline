#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * eDNA Full Ecological Analysis Pipeline — Standalone Entry Point
 *
 * Run this AFTER the main eDNA pipeline has completed, pointing
 * --results_dir at the main pipeline's output directory.
 *
 * Alternatively, supply ASV table + taxonomy directly per marker.
 *
 * Usage examples:
 *
 *   # From previous main pipeline run:
 *   nextflow run ecology_pipeline.nf \
 *       --results_dir results/ \
 *       --markers 16S,ITS \
 *       --metadata metadata.tsv \
 *       --outdir full_ecology/ \
 *       -profile singularity
 *
 *   # Direct file input (single marker):
 *   nextflow run ecology_pipeline.nf \
 *       --asv_table results/asv_tables/16S/16S.merged_asv_table.tsv \
 *       --taxonomy  results/asv_taxonomy/16S/16S.taxonomy_by_sequence.tsv \
 *       --marker    16S \
 *       --metadata  metadata.tsv \
 *       --outdir    full_ecology/ \
 *       -profile docker
 */

include { FULL_ECOLOGICAL_ANALYSIS } from './subworkflows/local/full_ecological_analysis'

// ─── Parameters ─────────────────────────────────────────────────────────────

params {
    // Input mode A: point at previous pipeline results dir
    results_dir       = null
    markers           = '16S'        // Comma-separated if results_dir mode

    // Input mode B: direct file inputs (single marker)
    asv_table         = null
    taxonomy          = null
    marker            = null

    // Common
    metadata          = null
    outdir            = 'full_ecology_results'

    // Ecology options
    ecology_group_var    = null   // Primary grouping variable name in metadata
    ecology_min_prevalence = 0.3  // Min fraction of samples ASV must appear in for network
    ecology_cor_cutoff   = 0.6    // Spearman |r| threshold for network edges

    max_memory = '128.GB'
    max_cpus   = 32
    max_time   = '48.h'
}

// ─── Input channel construction ──────────────────────────────────────────────

def buildInputChannel() {

    // Mode B: direct files for a single marker
    if (params.asv_table && params.taxonomy && params.marker) {
        return Channel.of([
            params.marker.toUpperCase(),
            file(params.asv_table, checkIfExists: true),
            file(params.taxonomy,  checkIfExists: true)
        ])
    }

    // Mode A: scan previous pipeline results directory
    if (!params.results_dir) {
        error "Provide either --results_dir or (--asv_table + --taxonomy + --marker)"
    }

    def results_dir = file(params.results_dir)
    if (!results_dir.exists()) {
        error "results_dir does not exist: ${params.results_dir}"
    }

    def markers_list = params.markers instanceof List
        ? params.markers
        : params.markers.tokenize(',').collect { it.trim().toUpperCase() }

    // Expected structure from main pipeline:
    //   {results_dir}/asv_tables/{MARKER}/{MARKER}.merged_asv_table.tsv
    //   {results_dir}/asv_taxonomy/{MARKER}/{MARKER}.taxonomy_by_sequence.tsv
    //
    // Must be the sequence-keyed file, not taxonomy/{MARKER}/{MARKER}.taxonomy.tsv
    // (that one's keyed by ASV label, not sequence — every ecology module joins
    // on the raw ASV sequence to match merged_asv_table.tsv's asv_id).

    def inputs = markers_list.collect { marker ->
        def asv_pattern = "${results_dir}/asv_tables/${marker}/**merged_asv_table.tsv"
        def tax_pattern = "${results_dir}/asv_taxonomy/${marker}/**taxonomy_by_sequence.tsv"

        def asv_files = file(asv_pattern)
        def tax_files = file(tax_pattern)

        // Handle glob returning list vs single file
        def asv_file = asv_files instanceof List
            ? (asv_files.isEmpty() ? null : asv_files[0])
            : (asv_files.exists() ? asv_files : null)

        def tax_file = tax_files instanceof List
            ? (tax_files.isEmpty() ? null : tax_files[0])
            : (tax_files.exists() ? tax_files : null)

        if (!asv_file) {
            log.warn "No merged ASV table found for ${marker} in ${results_dir}/asv_tables/${marker}/ — skipping."
            return null
        }
        if (!tax_file) {
            log.warn "No taxonomy table found for ${marker} in ${results_dir}/asv_taxonomy/${marker}/ — skipping."
            return null
        }

        log.info "Found ${marker}: ASV table=${asv_file.name}, taxonomy=${tax_file.name}"
        [ marker, asv_file, tax_file ]
    }
    .findAll { it != null }

    if (inputs.isEmpty()) {
        error "No valid marker data found in ${params.results_dir}. Check --markers and directory structure."
    }

    return Channel.fromList(inputs)
}

// ─── Workflow ────────────────────────────────────────────────────────────────

workflow {

    ch_ecology_input = buildInputChannel()

    def meta_file = params.metadata ? file(params.metadata, checkIfExists: true) : []

    FULL_ECOLOGICAL_ANALYSIS(
        ch_ecology_input,
        meta_file
    )
}

workflow.onComplete {
    log.info ""
    log.info "Full Ecological Analysis completed!"
    log.info "Duration : ${workflow.duration}"
    log.info "Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Results  : ${params.outdir}"
    log.info ""
    if (workflow.success) {
        log.info "Output structure:"
        log.info "  full_ecology/{MARKER}/00_report/           — HTML summary report"
        log.info "  full_ecology/{MARKER}/01_alpha_diversity/  — richness, diversity, rarefaction"
        log.info "  full_ecology/{MARKER}/02_beta_diversity/   — PERMANOVA, ANOSIM, distances"
        log.info "  full_ecology/{MARKER}/03_ordination/       — PCoA, NMDS, PCA, RDA, CCA"
        log.info "  full_ecology/{MARKER}/04_differential/     — DESeq2, ALDEx2"
        log.info "  full_ecology/{MARKER}/05_co_occurrence/    — network, hubs"
        log.info "  full_ecology/{MARKER}/06_indicator/        — IndVal, SIMPER, core"
        log.info "  full_ecology/{MARKER}/07_env_drivers/      — envfit, Mantel, varpart"
        log.info "  full_ecology/cross_marker/                 — inter-marker comparisons"
        log.info ""
    }
}
