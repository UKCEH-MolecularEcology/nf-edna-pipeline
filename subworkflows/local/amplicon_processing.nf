include { CUTADAPT           } from '../../modules/local/cutadapt/main'
include { DADA2_LEARN_ERRORS } from '../../modules/local/dada2_learn_errors/main'
include { DADA2_DENOISE      } from '../../modules/local/dada2_denoise/main'
include { VSEARCH_CHIMERA    } from '../../modules/local/vsearch_chimera/main'
include { TAXONOMY           } from '../../modules/local/taxonomy/main'

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

    // 2. Learn DADA2 error models per run (group samples by sequencing run)
    // If no run_id in meta, treat all samples as one run
    ch_trimmed_by_run = CUTADAPT.out.reads
        .map { meta, reads ->
            def run_id = meta.run ?: 'run1'
            [ run_id, meta, reads instanceof List ? reads : [reads] ]
        }
        .groupTuple(by: 0)
        .map { run_id, metas, reads_list ->
            // flatten [[R1,R2],[R1,R2],...] → [R1,R2,R1,R2,...] for path staging
            [ run_id, metas, reads_list.flatten() ]
        }

    DADA2_LEARN_ERRORS(
        ch_trimmed_by_run,
        marker,
        marker_params.trunc_len_f,
        marker_params.trunc_len_r,
        marker_params.max_ee_f,
        marker_params.max_ee_r
    )
    ch_versions = ch_versions.mix(DADA2_LEARN_ERRORS.out.versions.first())

    // 3. Denoise: sample inference + merge paired reads + make ASV table
    ch_denoise_input = CUTADAPT.out.reads
        .map { meta, reads ->
            def run_id = meta.run ?: 'run1'
            [ run_id, meta, reads ]
        }
        .combine(DADA2_LEARN_ERRORS.out.error_model, by: 0)

    DADA2_DENOISE(
        ch_denoise_input,
        marker,
        marker_params.trunc_len_f,
        marker_params.trunc_len_r,
        marker_params.max_ee_f,
        marker_params.max_ee_r,
        params.dada2_pool,
        marker_params.min_length,
        marker_params.max_length
    )
    ch_versions = ch_versions.mix(DADA2_DENOISE.out.versions.first())

    // 4. Chimera detection and removal (VSEARCH de novo + reference)
    VSEARCH_CHIMERA(
        DADA2_DENOISE.out.asv_seqs,
        marker
    )
    ch_versions = ch_versions.mix(VSEARCH_CHIMERA.out.versions.first())

    // 5. Filter ASV table to chimera-free ASVs and assign taxonomy
    TAXONOMY(
        VSEARCH_CHIMERA.out.nonchimeras,
        marker_params.tax_db,
        marker_params.tax_db_type,
        marker_params.tax_method,
        marker
    )
    ch_versions = ch_versions.mix(TAXONOMY.out.versions.first())

    emit:
    asv_table    = DADA2_DENOISE.out.asv_table      // [ meta, asv_table.rds ]
    asv_seqs     = VSEARCH_CHIMERA.out.nonchimeras  // [ meta, asv_seqs.fasta ]
    taxonomy     = TAXONOMY.out.taxonomy            // [ meta, taxonomy.tsv ]
    cutadapt_log = CUTADAPT.out.log
    versions     = ch_versions
}
