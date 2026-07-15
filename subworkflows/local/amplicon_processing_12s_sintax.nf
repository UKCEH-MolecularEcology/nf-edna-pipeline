include { SINTAX_DB_PREP         } from '../../modules/local/sintax_db_prep/main'
include { SINTAX_RUN             } from '../../modules/local/sintax_run/main'
include { SINTAX_MERGE           } from '../../modules/local/sintax_merge/main'
include { SINTAX_CLARE_ABUNDANCE } from '../../modules/local/sintax_clare_abundance/main'
include { SINTAX_BLANK_CLEANUP   } from '../../modules/local/sintax_blank_cleanup/main'

// 12S SINTAX taxonomy (against 3 reference databases) + LOD blank cleanup,
// ported from 12S-edna-dada2-tapirs-workflow. Runs on the marker-level
// MERGED, chimera-free ASV set for 12S (i.e. after MERGE_ASV_TABLES), not
// per-sample — this is the only place a sequence-keyed ASV set exists in
// this pipeline, matching what the Snakemake pipeline's SINTAX step expects.
workflow SINTAX_12S {

    take:
    ch_merged   // [ marker, merged_table, merged_fasta, asv_lookup ] — 12S only, singleton

    main:
    ch_versions = Channel.empty()
    s = params.sintax

    def missing = s.databases.findAll { !s.db_paths[it] || !file(s.db_paths[it]).exists() }
    if (missing) {
        error "params.sintax.enabled is true but reference database path(s) missing/unset for: ${missing.join(', ')}"
    }

    // Normalize each reference DB to a plain FASTA (gunzip MIDORI etc.) —
    // runs once per db; Nextflow's own -resume caching replaces the
    // Snakemake pipeline's manual unzip-once caching hack for free.
    ch_dbs = Channel.fromList(s.databases).map { db -> [db, file(s.db_paths[db])] }
    SINTAX_DB_PREP(ch_dbs)
    ch_ref_fastas = SINTAX_DB_PREP.out.fasta

    ch_combined = ch_merged
        .map { marker, table, fasta, lookup -> [marker, fasta, lookup, table] }
        .combine(ch_ref_fastas)
        .multiMap { marker, asv_fasta, asv_lookup, merged_table, db_name, ref_fasta ->
            merged_arg: [marker, asv_fasta, asv_lookup, merged_table]
            db_arg:     [db_name, ref_fasta]
        }

    SINTAX_RUN(ch_combined.merged_arg, ch_combined.db_arg, s.cutoff)
    ch_versions = ch_versions.mix(SINTAX_RUN.out.versions.first())

    // 3-way join across databases → comparison + summary tables
    ch_parsed_grouped = SINTAX_RUN.out.parsed
        .map { marker, db_name, tsv -> [marker, db_name, tsv] }
        .groupTuple(by: 0)

    ch_merge_input = ch_merged
        .map { marker, table, fasta, lookup -> [marker, lookup] }
        .join(ch_parsed_grouped, by: 0)

    SINTAX_MERGE(ch_merge_input)
    ch_versions = ch_versions.mix(SINTAX_MERGE.out.versions)

    // CLARE-specific abundance table feeds blank cleanup
    ch_clare_abundance = SINTAX_RUN.out.taxonomy_abundance
        .filter { marker, db_name, csv -> db_name == 'CLARE' }
        .map { marker, db_name, csv -> [marker, csv] }

    SINTAX_CLARE_ABUNDANCE(ch_clare_abundance)
    ch_versions = ch_versions.mix(SINTAX_CLARE_ABUNDANCE.out.versions)

    ch_blank_input = ch_clare_abundance.join(SINTAX_CLARE_ABUNDANCE.out.ncl_matrix, by: 0)
    SINTAX_BLANK_CLEANUP(ch_blank_input)
    ch_versions = ch_versions.mix(SINTAX_BLANK_CLEANUP.out.versions)

    emit:
    compare            = SINTAX_MERGE.out.compare
    database_summary   = SINTAX_MERGE.out.summary
    blank_cleaned_wide = SINTAX_BLANK_CLEANUP.out.cleaned_wide   // [ marker, ncl_cleaned_bothLOD.csv ]
    blank_cleaned_long = SINTAX_BLANK_CLEANUP.out.cleaned_long
    versions           = ch_versions
}
