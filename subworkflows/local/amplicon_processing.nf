include { CUTADAPT           } from '../../modules/local/cutadapt/main'
include { DADA2_FILTER       } from '../../modules/local/dada2_filter/main'
include { DADA2_LEARN_ERRORS } from '../../modules/local/dada2_learn_errors/main'
include { DADA2_DENOISE      } from '../../modules/local/dada2_denoise/main'
include { VSEARCH_CHIMERA    } from '../../modules/local/vsearch_chimera/main'

workflow AMPLICON_PROCESSING {

    take:
    ch_reads        // [ meta, reads ]
    marker          // String: '16S' | '18S' | 'ITS' | 'CO1' | '12S'
    marker_params   // Map of primer/DB/QC params

    main:
    ch_versions = Channel.empty()

    // 1. Primer trimming with marker-specific primers
    ch_reads_with_primers = ch_reads.map { meta, reads ->
        def new_meta = meta + [
            fwd_primer:  marker_params.fwd_primer,
            rev_primer:  marker_params.rev_primer,
            min_length:  marker_params.min_length,
            max_length:  marker_params.max_length
        ]
        [ new_meta, reads ]
    }
    CUTADAPT(ch_reads_with_primers)
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    // 2. Filter each sample independently, as its own Nextflow task --
    // real parallelism across samples (separate OS processes), unlike
    // relying on filterAndTrim's own multithread param, which parallelises
    // by forking (mclapply) across the files handed to ONE call. Forking
    // has nothing to gain once filtering already happens one sample at a
    // time, and (separately) fork-based parallelism was confirmed to give
    // no speedup at all in this container.
    ch_filter_input = CUTADAPT.out.reads.map { meta, reads ->
        def run_id = meta.run ?: 'run1'
        [ run_id, meta, reads ]
    }

    DADA2_FILTER(
        ch_filter_input,
        marker,
        marker_params.trunc_len_f,
        marker_params.trunc_len_r,
        marker_params.max_ee_f,
        marker_params.max_ee_r
    )
    ch_versions = ch_versions.mix(DADA2_FILTER.out.versions.first())

    // 3. Learn DADA2 error models per run, from the already-filtered reads
    // (group by run so every sample in a sequencing run contributes to one
    // shared error model, same as before).
    ch_filtered_by_run = DADA2_FILTER.out.filtered
        .map { run_id, meta, filtered -> [ run_id, meta, filtered instanceof List ? filtered : [filtered] ] }
        .groupTuple(by: 0)
        .map { run_id, metas, filtered_list -> [ run_id, metas, filtered_list.flatten() ] }

    DADA2_LEARN_ERRORS(
        ch_filtered_by_run,
        marker
    )
    ch_versions = ch_versions.mix(DADA2_LEARN_ERRORS.out.versions.first())

    // 4. Denoise: sample inference + merge paired reads + make ASV table.
    // Reuses each sample's own DADA2_FILTER output directly -- no
    // re-filtering here.
    ch_all_filtered_by_run = DADA2_FILTER.out.filtered
        .map { run_id, meta, filtered -> [ run_id, filtered instanceof List ? filtered : [filtered] ] }
        .groupTuple(by: 0)
        .map { run_id, filtered_list -> [ run_id, filtered_list.flatten() ] }

    ch_all_stats_by_run = DADA2_FILTER.out.stats
        .map { run_id, meta, stats -> [ run_id, stats instanceof List ? stats : [stats] ] }
        .groupTuple(by: 0)
        .map { run_id, stats_list -> [ run_id, stats_list.flatten() ] }

    ch_denoise_input = DADA2_FILTER.out.filtered
        .map { run_id, meta, filtered -> [ run_id, meta ] }
        .combine(DADA2_LEARN_ERRORS.out.error_model, by: 0)
        .combine(ch_all_filtered_by_run,             by: 0)
        .combine(ch_all_stats_by_run,                by: 0)

    DADA2_DENOISE(
        ch_denoise_input,
        marker,
        params.dada2_pool,
        marker_params.min_length,
        marker_params.max_length
    )
    ch_versions = ch_versions.mix(DADA2_DENOISE.out.versions.first())

    // 5. Chimera detection and removal (VSEARCH de novo + reference)
    VSEARCH_CHIMERA(
        DADA2_DENOISE.out.asv_seqs,
        marker
    )
    ch_versions = ch_versions.mix(VSEARCH_CHIMERA.out.versions.first())

    // Taxonomy is no longer assigned here (per-sample, on private pre-merge
    // ASVs) — it now runs once, collectively, on the marker-level merged
    // ASV set after MERGE_ASV_TABLES. See subworkflows/local/collective_taxonomy.nf
    // (and amplicon_processing_12s_sintax.nf for 12S specifically).

    emit:
    asv_table    = DADA2_DENOISE.out.asv_table      // [ meta, asv_table.rds ]
    asv_seqs     = VSEARCH_CHIMERA.out.nonchimeras  // [ meta, asv_seqs.fasta ]
    cutadapt_log = CUTADAPT.out.log
    versions     = ch_versions
}
