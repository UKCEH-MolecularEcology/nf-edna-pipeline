#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * eDNA Metabarcoding Pipeline
 * Supports: 16S, 18S, ITS, CO1, 12S
 * Raw reads → QC → Primer trimming → DADA2 denoising → Taxonomy → Ecology
 */

include { AMPLICON_QC              } from './subworkflows/local/amplicon_qc'
include { AMPLICON_PROCESSING      } from './subworkflows/local/amplicon_processing'
include { ECOLOGICAL_ANALYSIS      } from './subworkflows/local/ecological_analysis'
include { FULL_ECOLOGICAL_ANALYSIS } from './subworkflows/local/full_ecological_analysis'
include { MERGE_ASV_TABLES         } from './modules/local/merge_asvtables/main'
include { MULTIQC                  } from './modules/local/multiqc/main'

// ─── Validate parameters ────────────────────────────────────────────────────

def valid_markers = ['16S', '18S', 'ITS', 'CO1', '12S']

def checkParams() {
    if (!params.input) {
        error "Please provide a samplesheet with --input"
    }
    if (!params.outdir) {
        error "Please provide an output directory with --outdir"
    }
    def markers = params.markers instanceof List ? params.markers : params.markers.tokenize(',')
    markers.each { m ->
        if (!valid_markers.contains(m.trim())) {
            error "Invalid marker '${m}'. Valid markers: ${valid_markers.join(', ')}"
        }
    }
}

// ─── Parse samplesheet ──────────────────────────────────────────────────────

def parseSamplesheet(csv) {
    Channel.fromPath(csv)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:     row.sample,
                marker: row.marker.toUpperCase(),
                single_end: row.fastq_2 ? false : true
            ]
            def reads = row.fastq_2
                ? [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true) ]
            [ meta, reads ]
        }
}

// ─── Main workflow ───────────────────────────────────────────────────────────

workflow {

    checkParams()

    def markers_list = params.markers instanceof List
        ? params.markers
        : params.markers.tokenize(',').collect { it.trim() }

    // Parse samplesheet → channel of [meta, reads]
    ch_reads = parseSamplesheet(params.input)

    // Filter to requested markers
    ch_reads_filtered = ch_reads
        .filter { meta, reads -> markers_list.contains(meta.marker) }

    // ── QC ────────────────────────────────────────────────────────────────
    AMPLICON_QC(ch_reads_filtered)

    // ── Per-marker processing ─────────────────────────────────────────────
    // Branch by marker so each gets its own primer params and DB
    ch_by_marker = ch_reads_filtered
        .branch {
            s16:  it[0].marker == '16S'
            s18:  it[0].marker == '18S'
            its:  it[0].marker == 'ITS'
            co1:  it[0].marker == 'CO1'
            s12:  it[0].marker == '12S'
        }

    ch_asv_tables    = Channel.empty()
    ch_asv_seqs      = Channel.empty()
    ch_taxonomy_tbls = Channel.empty()

    markers_list.each { marker ->
        def ch_marker_reads = marker == '16S'  ? ch_by_marker.s16
                            : marker == '18S'  ? ch_by_marker.s18
                            : marker == 'ITS'  ? ch_by_marker.its
                            : marker == 'CO1'  ? ch_by_marker.co1
                            : ch_by_marker.s12

        def marker_params = loadMarkerParams(marker)

        AMPLICON_PROCESSING(
            ch_marker_reads,
            marker,
            marker_params
        )

        ch_asv_tables    = ch_asv_tables.mix(AMPLICON_PROCESSING.out.asv_table)
        ch_asv_seqs      = ch_asv_seqs.mix(AMPLICON_PROCESSING.out.asv_seqs)
        ch_taxonomy_tbls = ch_taxonomy_tbls.mix(AMPLICON_PROCESSING.out.taxonomy)
    }

    // ── Merge per-sample ASV tables → per-marker combined table ───────────
    ch_asv_grouped = ch_asv_tables
        .map { meta, tbl -> [ meta.marker, meta, tbl ] }
        .groupTuple(by: 0)

    MERGE_ASV_TABLES(ch_asv_grouped)

    // ── MultiQC report ────────────────────────────────────────────────────
    ch_multiqc_files = AMPLICON_QC.out.fastqc_zip
        .mix(AMPLICON_QC.out.cutadapt_log)
        .collect()
    MULTIQC(ch_multiqc_files)

    // ── Ecological analysis ───────────────────────────────────────────────
    // Build shared input channel: [ marker, asv_table, taxonomy ]
    ch_taxonomy_by_marker = ch_taxonomy_tbls
        .map { meta, tbl -> [ meta.marker, tbl ] }
        .groupTuple()
        .map { marker, tbls -> [ marker, tbls[0] ] }

    ch_ecology_input = MERGE_ASV_TABLES.out.merged_table
        .join(ch_taxonomy_by_marker, by: 0)

    def meta_file = params.metadata ? file(params.metadata) : []

    // Basic ecology (barplots, basic PCoA, diversity)
    if (params.run_ecology) {
        ECOLOGICAL_ANALYSIS(ch_ecology_input, meta_file)
    }

    // Full ecological analysis suite (runs after basic ecology)
    if (params.run_full_ecology) {
        FULL_ECOLOGICAL_ANALYSIS(ch_ecology_input, meta_file)
    }
}

// ─── Helper: load marker-specific primer and DB params ──────────────────────

def loadMarkerParams(marker) {
    def primers = params.primers[marker]
    def db      = params.databases[marker]

    if (!primers) error "No primer configuration found for marker: ${marker}"
    if (!db)      error "No database configuration found for marker: ${marker}"

    return [
        fwd_primer:    primers.fwd,
        rev_primer:    primers.rev,
        min_length:    primers.min_length ?: 50,
        max_length:    primers.max_length ?: 600,
        trunc_len_f:   primers.trunc_len_f ?: 0,
        trunc_len_r:   primers.trunc_len_r ?: 0,
        max_ee_f:      primers.max_ee_f ?: 2,
        max_ee_r:      primers.max_ee_r ?: 2,
        tax_db:        file(db.path, checkIfExists: true),
        tax_db_type:   db.type,   // silva | unite | midori | pr2
        tax_method:    db.method  // dada2 | blast | kraken2
    ]
}

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed!"
    log.info "Duration : ${workflow.duration}"
    log.info "Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Results  : ${params.outdir}"
    log.info ""
}
