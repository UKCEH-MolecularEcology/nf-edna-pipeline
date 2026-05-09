#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * eDNA Metabarcoding Pipeline
 * Supports: 16S, 18S, ITS, CO1, 12S, RBCL
 * Raw reads → QC → Primer trimming → DADA2 denoising → Taxonomy → Ecology
 *
 * Input options (mutually exclusive):
 *   --input      Samplesheet CSV (sample, fastq_1, fastq_2, marker)
 *   --fastq_dir  Directory of FASTQ files; markers auto-detected from filenames
 */

include { AMPLICON_QC              } from './subworkflows/local/amplicon_qc'
include { AMPLICON_PROCESSING      } from './subworkflows/local/amplicon_processing'
include { ECOLOGICAL_ANALYSIS      } from './subworkflows/local/ecological_analysis'
include { FULL_ECOLOGICAL_ANALYSIS } from './subworkflows/local/full_ecological_analysis'
include { MERGE_ASV_TABLES         } from './modules/local/merge_asvtables/main'
include { MULTIQC                  } from './modules/local/multiqc/main'

// ─── Validate parameters ────────────────────────────────────────────────────

def checkParams() {
    // Defined inside the function so it is in scope for DSL2 method compilation
    def valid_markers = ['16S', '18S', 'ITS', 'CO1', '12S', 'RBCL']

    if (!params.input && !params.fastq_dir) {
        error "Provide either --input (samplesheet CSV) or --fastq_dir (FASTQ directory)"
    }
    if (params.input && params.fastq_dir) {
        error "Specify either --input or --fastq_dir, not both"
    }
    if (!params.outdir) {
        error "Please provide an output directory with --outdir"
    }
    // Marker validation only applies when using --input (fastq_dir auto-detects them)
    if (params.input) {
        def markers = params.markers instanceof List ? params.markers : params.markers.tokenize(',')
        markers.each { m ->
            if (!valid_markers.contains(m.trim())) {
                error "Invalid marker '${m}'. Valid markers: ${valid_markers.join(', ')}"
            }
        }
    }
}

// ─── Parse samplesheet ──────────────────────────────────────────────────────

def parseSamplesheet(csv) {
    Channel.fromPath(csv)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:         row.sample,
                marker:     row.marker.toUpperCase(),
                single_end: row.fastq_2 ? false : true
            ]
            def reads = row.fastq_2
                ? [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true) ]
            [ meta, reads ]
        }
}

// ─── Auto-detect samples from FASTQ directory ───────────────────────────────

def parseFastqDir(fastq_dir) {
    // Defined inside the function so it is in scope for DSL2 method compilation
    def marker_aliases = [
        '16S': '16S', '18S': '18S', 'ITS': 'ITS', 'ITS1': 'ITS', 'ITS2': 'ITS',
        'CO1': 'CO1', 'COI': 'CO1', '12S': '12S', 'RBCL': 'RBCL'
    ]

    def dir = file(fastq_dir)
    if (!dir.isDirectory()) error "--fastq_dir is not a directory: ${fastq_dir}"

    def samples  = [:]   // [sample_id, marker] → [R1: file, R2: file]
    def skipped  = []

    dir.listFiles()
        .findAll { it.name.endsWith('.fastq.gz') }
        .sort    { it.name }
        .each    { f ->
            def m = (f.name =~ /_(S\d+)_(R[12])_001\.fastq\.gz$/)
            if (!m) { skipped << f.name; return }

            def read = m[0][2].toUpperCase()
            def stem = f.name[0..(f.name.size() - m[0][0].size() - 1)]

            def sample_id = null
            def marker    = null
            for (sep in ['-', '_']) {
                def idx = stem.lastIndexOf(sep as String)
                if (idx < 0) continue
                def candidate = stem[(idx + 1)..-1].toUpperCase()
                if (marker_aliases.containsKey(candidate)) {
                    sample_id = stem[0..(idx - 1)]
                    marker    = marker_aliases[candidate]
                    break
                }
            }
            if (!sample_id) { skipped << f.name; return }

            def key = [sample_id, marker]
            if (!samples.containsKey(key)) samples[key] = [:]
            if (samples[key].containsKey(read)) {
                error "Duplicate ${read} for sample='${sample_id}' marker='${marker}': ${f}"
            }
            samples[key][read] = f
        }

    if (skipped) {
        log.warn "Skipped ${skipped.size()} unrecognised FASTQ file(s):\n  ${skipped.join('\n  ')}"
    }
    if (!samples) error "No samples detected in: ${fastq_dir}"
    return samples
}

// ─── Write auto-generated samplesheet ───────────────────────────────────────

def writeSamplesheet(samples, path_str) {
    def out = file(path_str)
    out.parent.mkdirs()
    def lines = ['sample,fastq_1,fastq_2,marker']
    samples.sort { a, b -> (a.key[0] <=> b.key[0]) ?: (a.key[1] <=> b.key[1]) }
           .each { key, reads ->
                lines << "${key[0]},${reads.R1 ?: ''},${reads.R2 ?: ''},${key[1]}"
           }
    out.text = lines.join('\n') + '\n'
    log.info "Auto-generated samplesheet: ${out}"
}

// ─── Main workflow ───────────────────────────────────────────────────────────

workflow {

    checkParams()

    // ── Resolve input: samplesheet OR FASTQ directory ────────────────────
    def markers_list
    def ch_reads

    if (params.fastq_dir) {
        def detected = parseFastqDir(params.fastq_dir)
        markers_list = detected.keySet().collect { it[1] }.unique().sort()
        log.info "Auto-detected markers: ${markers_list.join(', ')}"
        writeSamplesheet(detected, "${params.outdir}/samplesheet_detected.csv")

        ch_reads = Channel.fromList(
            detected.collect { key, reads ->
                def meta = [
                    id:         key[0],
                    marker:     key[1],
                    single_end: !reads.containsKey('R2')
                ]
                def files = reads.containsKey('R2')
                    ? [ reads.R1, reads.R2 ]
                    : [ reads.R1 ]
                [ meta, files ]
            }
        )
    } else {
        markers_list = params.markers instanceof List
            ? params.markers
            : params.markers.tokenize(',').collect { it.trim() }
        ch_reads = parseSamplesheet(params.input)
    }

    // ── Metadata (NO_FILE placeholder when absent) ────────────────────────
    def meta_file = params.metadata
        ? file(params.metadata)
        : file("${projectDir}/assets/NO_FILE")

    // ── Filter to requested markers ───────────────────────────────────────
    ch_reads_filtered = ch_reads
        .filter { meta, reads -> markers_list.contains(meta.marker) }

    // ── QC ────────────────────────────────────────────────────────────────
    AMPLICON_QC(ch_reads_filtered)

    // ── Per-marker processing ─────────────────────────────────────────────
    ch_by_marker = ch_reads_filtered
        .branch {
            s16:  it[0].marker == '16S'
            s18:  it[0].marker == '18S'
            its:  it[0].marker == 'ITS'
            co1:  it[0].marker == 'CO1'
            s12:  it[0].marker == '12S'
            rbcl: it[0].marker == 'RBCL'
        }

    ch_asv_tables    = Channel.empty()
    ch_asv_seqs      = Channel.empty()
    ch_taxonomy_tbls = Channel.empty()

    markers_list.each { marker ->
        def ch_marker_reads = marker == '16S'  ? ch_by_marker.s16
                            : marker == '18S'  ? ch_by_marker.s18
                            : marker == 'ITS'  ? ch_by_marker.its
                            : marker == 'CO1'  ? ch_by_marker.co1
                            : marker == '12S'  ? ch_by_marker.s12
                            : ch_by_marker.rbcl

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
    ch_taxonomy_by_marker = ch_taxonomy_tbls
        .map { meta, tbl -> [ meta.marker, tbl ] }
        .groupTuple()
        .map { marker, tbls -> [ marker, tbls[0] ] }

    ch_ecology_input = MERGE_ASV_TABLES.out.merged_table
        .join(ch_taxonomy_by_marker, by: 0)

    if (params.run_ecology) {
        ECOLOGICAL_ANALYSIS(ch_ecology_input, meta_file)
    }

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

    def tax_db_file = file(db.path)
    if (!tax_db_file.exists()) {
        error """\
            Database not found for ${marker}: ${db.path}
            Run the download script first:
              bash ${projectDir}/assets/download_databases.sh ${projectDir}/databases/
            ITS (UNITE) requires a manual download — see assets/download_databases.sh for details.
            """.stripIndent()
    }

    return [
        fwd_primer:    primers.fwd,
        rev_primer:    primers.rev,
        min_length:    primers.min_length ?: 50,
        max_length:    primers.max_length ?: 600,
        trunc_len_f:   primers.trunc_len_f ?: 0,
        trunc_len_r:   primers.trunc_len_r ?: 0,
        max_ee_f:      primers.max_ee_f ?: 2,
        max_ee_r:      primers.max_ee_r ?: 2,
        tax_db:        tax_db_file,
        tax_db_type:   db.type,
        tax_method:    db.method
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
