include { HONEYPI_DOWNLOAD_DB      } from '../../modules/local/honeypi_download_db/main'
include { HONEYPI_TRIM_GALORE      } from '../../modules/local/honeypi_trim_galore/main'
include { HONEYPI_DADA2            } from '../../modules/local/honeypi_dada2/main'
include { HONEYPI_ITSX             } from '../../modules/local/honeypi_itsx/main'
include { HONEYPI_CONSOLIDATE      } from '../../modules/local/honeypi_consolidate/main'
include { HONEYPI_RDP_CLASSIFIER   } from '../../modules/local/honeypi_rdp_classifier/main'
include { HONEYPI_FILTER_ASV_TABLE } from '../../modules/local/honeypi_filter_asv_table/main'
include { HONEYPI_MERGE_DUPLICATES } from '../../modules/local/honeypi_merge_duplicates/main'
include { HONEYPI_SUMMARY          } from '../../modules/local/honeypi_summary/main'
include { HONEYPI_STANDARDIZE      } from '../../modules/local/honeypi_standardize/main'

// Plant ITS (pITS) processing via the honeypi workflow (Trim Galore -> joint
// DADA2 -> ITSx -> consolidate -> RDP classify -> merge), ported in from
// https://github.com/UKCEH-MolecularEcology/nf-honeypi as a subworkflow
// rather than re-implemented against this pipeline's per-sample/per-run
// DADA2 modules -- honeypi denoises all samples jointly in one step, a
// structurally different shape from AMPLICON_PROCESSING's per-run channel
// design, so wrapping its own modules is far lower-risk than rewriting them.
//
// honeypi requires hyphen-only sample IDs (its RDP/DADA2 steps reject
// underscores); this pipeline's sample IDs are underscore-based. IDs are
// sanitised on the way in and restored on the way out via a generated
// sample_id_map.tsv, so callers never see the hyphenated form.
workflow PLANT_ITS_PROCESSING {

    take:
    ch_reads   // [ meta(id, marker, ...), [ r1, r2 ] ]
    marker     // String, e.g. 'PITS'

    main:
    ch_versions = Channel.empty()

    ch_reads_sanitized = ch_reads.map { meta, reads ->
        def sano = (meta.id as String).replaceAll('_', '-').replaceAll(/-+$/, '')
        // Mirrors R's make.names(), which DADA2/data.frame apply to sample
        // (column) names: invalid leading char gets an 'X' prefix, hyphens
        // become dots. Precomputing this lets HONEYPI_STANDARDIZE translate
        // the R-mangled column names straight back to the original sample ID.
        def r_col = sano.replaceAll('-', '.')
        if (r_col =~ /^[0-9]/) r_col = 'X' + r_col
        def new_meta = meta + [orig_id: meta.id, id: sano, r_colname: r_col]
        [ new_meta, reads[0], reads[1] ]
    }

    ch_id_map = ch_reads_sanitized
        .map { meta, r1, r2 -> "${meta.r_colname}\t${meta.orig_id}" }
        .collectFile(name: 'sample_id_map.tsv', newLine: true, seed: 'r_colname\torig_id', storeDir: "${params.outdir}/asv_taxonomy/${marker}")

    // ── Per-sample trimming ───────────────────────────────────────────────
    HONEYPI_TRIM_GALORE(ch_reads_sanitized)
    ch_versions = ch_versions.mix(HONEYPI_TRIM_GALORE.out.versions.first())

    // ── Joint DADA2 denoising (all samples together, honeypi's own design) ─
    ch_r1_all = HONEYPI_TRIM_GALORE.out.reads.map { meta, r1, r2 -> r1 }.collect()
    ch_r2_all = HONEYPI_TRIM_GALORE.out.reads.map { meta, r1, r2 -> r2 }.collect()

    HONEYPI_DADA2(ch_r1_all, ch_r2_all)
    ch_versions = ch_versions.mix(HONEYPI_DADA2.out.versions)

    // ── ITS extraction ────────────────────────────────────────────────────
    HONEYPI_ITSX(
        HONEYPI_DADA2.out.asvs,
        params.honeypi.its_region,
        params.honeypi.its_min_len,
        params.honeypi.its_max_len
    )
    ch_versions = ch_versions.mix(HONEYPI_ITSX.out.versions)

    // ── Consolidate (ITSxed + not-ITSxed) ────────────────────────────────
    HONEYPI_CONSOLIDATE(HONEYPI_DADA2.out.asvs, HONEYPI_ITSX.out.itsxed)
    ch_versions = ch_versions.mix(HONEYPI_CONSOLIDATE.out.versions)

    // ── RDP database ──────────────────────────────────────────────────────
    if (params.honeypi.rdp_db_dir && file(params.honeypi.rdp_db_dir).exists()) {
        ch_db_dir = Channel.value(file(params.honeypi.rdp_db_dir))
    } else {
        def db_url = (params.honeypi.its_region == 'ITS1') ? params.honeypi.rdp_db_its1_url : params.honeypi.rdp_db_its2_url
        HONEYPI_DOWNLOAD_DB(Channel.value(db_url), Channel.value(params.honeypi.its_region))
        ch_db_dir = HONEYPI_DOWNLOAD_DB.out.db_dir
        ch_versions = ch_versions.mix(HONEYPI_DOWNLOAD_DB.out.versions)
    }

    // ── RDP classification ────────────────────────────────────────────────
    HONEYPI_RDP_CLASSIFIER(
        HONEYPI_CONSOLIDATE.out.fasta,
        ch_db_dir,
        params.honeypi.rdp_confidence
    )
    ch_versions = ch_versions.mix(HONEYPI_RDP_CLASSIFIER.out.versions)

    // ── Filter count table to classified ASVs, then merge same-taxonomy ASVs
    HONEYPI_FILTER_ASV_TABLE(HONEYPI_DADA2.out.counts, HONEYPI_CONSOLIDATE.out.fasta)
    ch_versions = ch_versions.mix(HONEYPI_FILTER_ASV_TABLE.out.versions)

    HONEYPI_MERGE_DUPLICATES(HONEYPI_FILTER_ASV_TABLE.out.counts, HONEYPI_RDP_CLASSIFIER.out.taxonomy)
    ch_versions = ch_versions.mix(HONEYPI_MERGE_DUPLICATES.out.versions)

    // ── Native honeypi_output summary folder (unchanged deliverable) ───────
    ch_error_rates = HONEYPI_DADA2.out.error_plots.ifEmpty(file("${projectDir}/assets/NO_FILE"))
    HONEYPI_SUMMARY(
        HONEYPI_CONSOLIDATE.out.fasta,
        HONEYPI_MERGE_DUPLICATES.out.counts,
        HONEYPI_FILTER_ASV_TABLE.out.counts,
        HONEYPI_RDP_CLASSIFIER.out.taxonomy,
        HONEYPI_DADA2.out.stats,
        ch_error_rates
    )

    // ── Standardize to this pipeline's sequence-keyed table shape, with
    //    original (underscore) sample IDs restored, so it plugs into the
    //    same MERGE_ASV_TABLES / ecology join as every other marker ────────
    ch_standardize_input = Channel.value(marker)
        .combine(HONEYPI_FILTER_ASV_TABLE.out.counts)
        .combine(HONEYPI_CONSOLIDATE.out.fasta)
        .combine(HONEYPI_RDP_CLASSIFIER.out.taxonomy)
        .combine(ch_id_map)

    HONEYPI_STANDARDIZE(ch_standardize_input)
    ch_versions = ch_versions.mix(HONEYPI_STANDARDIZE.out.versions)

    emit:
    merged_table         = HONEYPI_STANDARDIZE.out.merged_table          // [ marker, tsv ]
    taxonomy_for_ecology = HONEYPI_STANDARDIZE.out.taxonomy_by_sequence  // [ marker, tsv ]
    honeypi_counts       = HONEYPI_MERGE_DUPLICATES.out.counts
    versions             = ch_versions
}
