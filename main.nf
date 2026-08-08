#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * eDNA Metabarcoding Pipeline
 * Supports: 16S, 18S, ITS, CO1, 12S, RBCL, PITS (plant ITS2, via the honeypi
 * subworkflow -- ported from https://github.com/UKCEH-MolecularEcology/nf-honeypi)
 * Raw reads → QC → Primer trimming → DADA2 denoising → Taxonomy → Ecology
 *
 * Input options (mutually exclusive):
 *   --input      Samplesheet CSV (sample, fastq_1, fastq_2, marker)
 *   --fastq_dir  Directory of FASTQ files; markers auto-detected from filenames
 *                ('fITS' filename tokens map to fungal ITS/UNITE, 'pITS' tokens
 *                route to the PITS/honeypi branch)
 */

include { AMPLICON_QC                            } from './subworkflows/local/amplicon_qc'
include { AMPLICON_PROCESSING as AMPLICON_PROCESSING_16S  } from './subworkflows/local/amplicon_processing'
include { AMPLICON_PROCESSING as AMPLICON_PROCESSING_18S  } from './subworkflows/local/amplicon_processing'
include { AMPLICON_PROCESSING as AMPLICON_PROCESSING_ITS  } from './subworkflows/local/amplicon_processing'
include { AMPLICON_PROCESSING as AMPLICON_PROCESSING_CO1  } from './subworkflows/local/amplicon_processing'
include { AMPLICON_PROCESSING as AMPLICON_PROCESSING_12S  } from './subworkflows/local/amplicon_processing'
include { AMPLICON_PROCESSING as AMPLICON_PROCESSING_RBCL } from './subworkflows/local/amplicon_processing'
include { TAPIRS_12S                            } from './subworkflows/local/amplicon_processing_12s_tapirs'
include { SINTAX_12S                            } from './subworkflows/local/amplicon_processing_12s_sintax'
include { TAPIRS_BLAST_LCA_CO1                  } from './subworkflows/local/tapirs_blast_lca_co1'
include { COLLECTIVE_TAXONOMY as COLLECTIVE_TAXONOMY_16S  } from './subworkflows/local/collective_taxonomy'
include { COLLECTIVE_TAXONOMY as COLLECTIVE_TAXONOMY_18S  } from './subworkflows/local/collective_taxonomy'
include { COLLECTIVE_TAXONOMY as COLLECTIVE_TAXONOMY_ITS  } from './subworkflows/local/collective_taxonomy'
include { COLLECTIVE_TAXONOMY as COLLECTIVE_TAXONOMY_CO1  } from './subworkflows/local/collective_taxonomy'
include { COLLECTIVE_TAXONOMY as COLLECTIVE_TAXONOMY_RBCL } from './subworkflows/local/collective_taxonomy'
include { PLANT_ITS_PROCESSING                  } from './subworkflows/local/plant_its_processing'
include { ECOLOGICAL_ANALYSIS                  } from './subworkflows/local/ecological_analysis'
include { FULL_ECOLOGICAL_ANALYSIS             } from './subworkflows/local/full_ecological_analysis'
include { MERGE_ASV_TABLES                     } from './modules/local/merge_asvtables/main'
include { MULTIQC                              } from './modules/local/multiqc/main'
include { SETUP_DATABASES                      } from './modules/local/setup_databases/main'

// ─── Validate parameters ────────────────────────────────────────────────────

def checkParams() {
    // Defined inside the function so it is in scope for DSL2 method compilation
    def valid_markers = ['16S', '18S', 'ITS', 'CO1', '12S', 'RBCL', 'PITS']

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
        'CO1': 'CO1', 'COI': 'CO1', '12S': '12S', 'RBCL': 'RBCL',
        // Fungal ITS (fITS) reuses the generic ITS/UNITE path; plant ITS
        // (pITS) routes to PITS, processed by the honeypi subworkflow
        // (different pipeline shape + a plant-specific RDP database).
        'FITS': 'ITS', 'PITS': 'PITS'
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

    // ── Database setup (auto-download if any are missing) ─────────────────
    def missing_db_markers = markers_list.findAll { !file(params.databases[it].path).exists() }

    if (missing_db_markers) {
        if (params.skip_db_download) {
            error """\
                Missing databases for: ${missing_db_markers.join(', ')}
                Run:  bash ${projectDir}/assets/download_databases.sh ${projectDir}/databases/ ${missing_db_markers.join(' ')}
                Or omit --skip_db_download to let the pipeline download them automatically.
                ITS (UNITE) always requires a manual download — see assets/download_databases.sh.
                """.stripIndent()
        }
        log.info "Auto-downloading missing databases for: ${missing_db_markers.join(', ')}"
        SETUP_DATABASES(Channel.value(missing_db_markers))
    }

    // QC runs immediately; processing is gated on database readiness
    AMPLICON_QC(ch_reads_filtered)

    ch_proc_input = missing_db_markers
        ? ch_reads_filtered.combine(SETUP_DATABASES.out.ready).map { meta, reads, _ -> [meta, reads] }
        : ch_reads_filtered

    // ── Per-marker processing ─────────────────────────────────────────────
    ch_by_marker = ch_proc_input
        .branch {
            s16:  it[0].marker == '16S'
            s18:  it[0].marker == '18S'
            its:  it[0].marker == 'ITS'
            co1:  it[0].marker == 'CO1'
            s12:  it[0].marker == '12S'
            rbcl: it[0].marker == 'RBCL'
            pits: it[0].marker == 'PITS'
        }

    ch_asv_tables    = Channel.empty()
    ch_asv_seqs      = Channel.empty()
    ch_cutadapt_logs = Channel.empty()

    // DSL2 requires each subworkflow to be called at most once; use aliased imports
    if (markers_list.contains('16S')) {
        AMPLICON_PROCESSING_16S(ch_by_marker.s16, '16S', loadMarkerParams('16S'))
        ch_asv_tables    = ch_asv_tables.mix(AMPLICON_PROCESSING_16S.out.asv_table)
        ch_asv_seqs      = ch_asv_seqs.mix(AMPLICON_PROCESSING_16S.out.asv_seqs)
        ch_cutadapt_logs = ch_cutadapt_logs.mix(AMPLICON_PROCESSING_16S.out.cutadapt_log)
    }
    if (markers_list.contains('18S')) {
        AMPLICON_PROCESSING_18S(ch_by_marker.s18, '18S', loadMarkerParams('18S'))
        ch_asv_tables    = ch_asv_tables.mix(AMPLICON_PROCESSING_18S.out.asv_table)
        ch_asv_seqs      = ch_asv_seqs.mix(AMPLICON_PROCESSING_18S.out.asv_seqs)
        ch_cutadapt_logs = ch_cutadapt_logs.mix(AMPLICON_PROCESSING_18S.out.cutadapt_log)
    }
    if (markers_list.contains('ITS')) {
        AMPLICON_PROCESSING_ITS(ch_by_marker.its, 'ITS', loadMarkerParams('ITS'))
        ch_asv_tables    = ch_asv_tables.mix(AMPLICON_PROCESSING_ITS.out.asv_table)
        ch_asv_seqs      = ch_asv_seqs.mix(AMPLICON_PROCESSING_ITS.out.asv_seqs)
        ch_cutadapt_logs = ch_cutadapt_logs.mix(AMPLICON_PROCESSING_ITS.out.cutadapt_log)
    }
    if (markers_list.contains('CO1')) {
        AMPLICON_PROCESSING_CO1(ch_by_marker.co1, 'CO1', loadMarkerParams('CO1'))
        ch_asv_tables    = ch_asv_tables.mix(AMPLICON_PROCESSING_CO1.out.asv_table)
        ch_asv_seqs      = ch_asv_seqs.mix(AMPLICON_PROCESSING_CO1.out.asv_seqs)
        ch_cutadapt_logs = ch_cutadapt_logs.mix(AMPLICON_PROCESSING_CO1.out.cutadapt_log)
    }
    if (markers_list.contains('12S')) {
        AMPLICON_PROCESSING_12S(ch_by_marker.s12, '12S', loadMarkerParams('12S'))
        ch_asv_tables    = ch_asv_tables.mix(AMPLICON_PROCESSING_12S.out.asv_table)
        ch_asv_seqs      = ch_asv_seqs.mix(AMPLICON_PROCESSING_12S.out.asv_seqs)
        ch_cutadapt_logs = ch_cutadapt_logs.mix(AMPLICON_PROCESSING_12S.out.cutadapt_log)

        // Standalone Tapirs branch (fastp/vsearch/BLAST/Kraken2/MLCA) — fed
        // raw reads directly, self-contained, does not feed ecology.
        if (params.tapirs.enabled) {
            TAPIRS_12S(ch_by_marker.s12)
        }
    }
    if (markers_list.contains('RBCL')) {
        AMPLICON_PROCESSING_RBCL(ch_by_marker.rbcl, 'RBCL', loadMarkerParams('RBCL'))
        ch_asv_tables    = ch_asv_tables.mix(AMPLICON_PROCESSING_RBCL.out.asv_table)
        ch_asv_seqs      = ch_asv_seqs.mix(AMPLICON_PROCESSING_RBCL.out.asv_seqs)
        ch_cutadapt_logs = ch_cutadapt_logs.mix(AMPLICON_PROCESSING_RBCL.out.cutadapt_log)
    }

    // PITS (plant ITS2): structurally its own pipeline (joint DADA2 +
    // ITSx + RDP, ported from nf-honeypi) -- does not feed the generic
    // per-sample AMPLICON_PROCESSING/MERGE_ASV_TABLES/COLLECTIVE_TAXONOMY
    // path above. It emits an already-merged, sequence-keyed table plus a
    // taxonomy-by-sequence table in the same shape those produce, joined
    // in below (ch_merged_tables_all / ch_taxonomy_by_marker).
    if (markers_list.contains('PITS')) {
        PLANT_ITS_PROCESSING(ch_by_marker.pits, 'PITS', loadMarkerParams('PITS'))
        ch_cutadapt_logs = ch_cutadapt_logs.mix(PLANT_ITS_PROCESSING.out.cutadapt_log)
    }

    // ── Merge per-sample ASV tables → per-marker combined table ───────────
    ch_asv_grouped = ch_asv_tables
        .map { meta, tbl -> [ meta.marker, meta, tbl ] }
        .groupTuple(by: 0)

    MERGE_ASV_TABLES(ch_asv_grouped)

    // ── Collective (post-merge) taxonomy ───────────────────────────────────
    // Taxonomy is assigned once per marker, on the marker-level merged,
    // chimera-free ASV set — not per-sample, on private pre-merge ASVs (the
    // old behavior). This matches the source Snakemake pipeline's own
    // design (DADA2 + SINTAX both run collectively across all samples).
    ch_merged_all = MERGE_ASV_TABLES.out.merged_table
        .join(MERGE_ASV_TABLES.out.merged_fasta, by: 0)
        .join(MERGE_ASV_TABLES.out.asv_lookup, by: 0)

    // 12S: SINTAX taxonomy (3 databases) + LOD blank cleanup.
    if (markers_list.contains('12S') && params.sintax.enabled) {
        ch_merged_12s = ch_merged_all.filter { marker, table, fasta, lookup -> marker == '12S' }
        SINTAX_12S(ch_merged_12s)
    }

    // 16S/18S/ITS/CO1/RBCL: generic dada2/rdp/blast taxonomy (per marker's
    // params.databases[marker].method), joined with per-sample counts into
    // one combined ASV + sequence + taxonomy-lineage table.
    if (markers_list.contains('16S')) {
        COLLECTIVE_TAXONOMY_16S(ch_merged_all.filter { m, t, f, l -> m == '16S' }, '16S', loadMarkerParams('16S'))
    }
    if (markers_list.contains('18S')) {
        COLLECTIVE_TAXONOMY_18S(ch_merged_all.filter { m, t, f, l -> m == '18S' }, '18S', loadMarkerParams('18S'))
    }
    if (markers_list.contains('ITS')) {
        COLLECTIVE_TAXONOMY_ITS(ch_merged_all.filter { m, t, f, l -> m == 'ITS' }, 'ITS', loadMarkerParams('ITS'))
    }
    if (markers_list.contains('CO1')) {
        COLLECTIVE_TAXONOMY_CO1(ch_merged_all.filter { m, t, f, l -> m == 'CO1' }, 'CO1', loadMarkerParams('CO1'))

        // CO1's own BLAST + majority-vote-LCA branch (via coidb) -- a second,
        // independent taxonomy call for comparison, same shape as 12S's
        // Tapirs branch: standalone, does not feed ecology.
        if (params.tapirs_blast_lca_co1.enabled) {
            TAPIRS_BLAST_LCA_CO1(ch_merged_all.filter { m, t, f, l -> m == 'CO1' }, params.tapirs_blast_lca_co1)
        }
    }
    if (markers_list.contains('RBCL')) {
        COLLECTIVE_TAXONOMY_RBCL(ch_merged_all.filter { m, t, f, l -> m == 'RBCL' }, 'RBCL', loadMarkerParams('RBCL'))
    }

    // ── MultiQC report ────────────────────────────────────────────────────
    ch_multiqc_files = AMPLICON_QC.out.fastqc_zip.map { meta, zip -> zip }
        .mix(ch_cutadapt_logs.map { meta, log -> log })
        .collect()
    MULTIQC(ch_multiqc_files)

    // ── Ecological analysis ───────────────────────────────────────────────
    // Every marker's taxonomy now comes from its own collective (post-merge)
    // step, keyed by sequence to match merged_asv_table.tsv's asv_id — no
    // more picking one arbitrary sample's taxonomy as a marker-wide stand-in.
    ch_taxonomy_by_marker = Channel.empty()
    if (markers_list.contains('12S') && params.sintax.enabled) {
        ch_taxonomy_by_marker = ch_taxonomy_by_marker.mix(SINTAX_12S.out.blank_cleaned_wide)
    }
    if (markers_list.contains('16S')) {
        ch_taxonomy_by_marker = ch_taxonomy_by_marker.mix(COLLECTIVE_TAXONOMY_16S.out.taxonomy_for_ecology)
    }
    if (markers_list.contains('18S')) {
        ch_taxonomy_by_marker = ch_taxonomy_by_marker.mix(COLLECTIVE_TAXONOMY_18S.out.taxonomy_for_ecology)
    }
    if (markers_list.contains('ITS')) {
        ch_taxonomy_by_marker = ch_taxonomy_by_marker.mix(COLLECTIVE_TAXONOMY_ITS.out.taxonomy_for_ecology)
    }
    if (markers_list.contains('CO1')) {
        ch_taxonomy_by_marker = ch_taxonomy_by_marker.mix(COLLECTIVE_TAXONOMY_CO1.out.taxonomy_for_ecology)
    }
    if (markers_list.contains('RBCL')) {
        ch_taxonomy_by_marker = ch_taxonomy_by_marker.mix(COLLECTIVE_TAXONOMY_RBCL.out.taxonomy_for_ecology)
    }

    // PITS bypasses MERGE_ASV_TABLES/COLLECTIVE_TAXONOMY entirely (see
    // dispatch above) but emits the same [ marker, tsv ] shapes directly,
    // so it joins the ecology input alongside every other marker's output.
    ch_merged_tables_all = MERGE_ASV_TABLES.out.merged_table
    if (markers_list.contains('PITS')) {
        ch_merged_tables_all = ch_merged_tables_all.mix(PLANT_ITS_PROCESSING.out.merged_table)
        ch_taxonomy_by_marker = ch_taxonomy_by_marker.mix(PLANT_ITS_PROCESSING.out.taxonomy_for_ecology)
    }

    ch_ecology_input = ch_merged_tables_all
        .join(ch_taxonomy_by_marker, by: 0)

    // Filter to markers whose ecology is not disabled via --skip_ecology_markers
    def skip_eco = params.skip_ecology_markers
        ? (params.skip_ecology_markers instanceof List
            ? params.skip_ecology_markers
            : params.skip_ecology_markers.tokenize(',').collect { it.trim().toUpperCase() })
        : []

    ch_ecology_filtered = skip_eco
        ? ch_ecology_input.filter { marker, tbl, tax -> !skip_eco.contains(marker.toUpperCase()) }
        : ch_ecology_input

    if (params.run_ecology) {
        ECOLOGICAL_ANALYSIS(ch_ecology_filtered, meta_file)
    }

    if (params.run_full_ecology) {
        FULL_ECOLOGICAL_ANALYSIS(ch_ecology_filtered, meta_file)
    }
}

// ─── Helper: load marker-specific primer and DB params ──────────────────────

def loadMarkerParams(marker) {
    def primers = params.primers[marker]
    def db      = params.databases[marker]

    if (!primers) error "No primer configuration found for marker: ${marker}"
    if (!db)      error "No database configuration found for marker: ${marker}"

    // 12S fusion-tag support: when a fusion tag was used in library prep,
    // the forward primer's first base appears as a degenerate R (A|G) in
    // the reads rather than the plain A.
    def fwd = primers.fwd
    if (marker == '12S' && primers.fusion_tag) {
        fwd = 'R' + fwd.substring(1)
    }

    return [
        fwd_primer:    fwd,
        rev_primer:    primers.rev,
        min_length:    primers.min_length ?: 50,
        max_length:    primers.max_length ?: 600,
        trunc_len_f:   primers.trunc_len_f ?: 0,
        trunc_len_r:   primers.trunc_len_r ?: 0,
        max_ee_f:      primers.max_ee_f ?: 2,
        max_ee_r:      primers.max_ee_r ?: 2,
        tax_db:        file(db.path),
        tax_db_type:   db.type,
        tax_method:    db.method,
        addspecies_db: db.addspecies_path ? file(db.addspecies_path) : file("${projectDir}/assets/NO_FILE")
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
