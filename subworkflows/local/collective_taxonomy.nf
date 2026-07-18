include { TAXONOMY           } from '../../modules/local/taxonomy/main'
include { ASV_TAXONOMY_TABLE } from '../../modules/local/asv_taxonomy_table/main'

// Collective (post-merge) taxonomy for a single marker: assigns taxonomy
// once to the marker-level merged, chimera-free ASV set (not per-sample,
// which is what used to happen here), then joins it with the merged
// abundance table into one combined output. Mirrors
// amplicon_processing_12s_sintax.nf's SINTAX_12S pattern, generalized to
// the dada2/rdp/blast methods used by the other 5 markers.
workflow COLLECTIVE_TAXONOMY {

    take:
    ch_merged       // [ marker, merged_table, merged_fasta, asv_lookup ] — singleton, one marker
    marker
    marker_params   // { tax_db, tax_db_type, tax_method, addspecies_db } from loadMarkerParams(marker)

    main:
    ch_versions = Channel.empty()

    TAXONOMY(
        ch_merged.map { m, table, fasta, lookup -> [ m, fasta ] },
        marker_params.tax_db,
        marker_params.tax_db_type,
        marker_params.tax_method,
        marker_params.addspecies_db
    )
    ch_versions = ch_versions.mix(TAXONOMY.out.versions)

    ch_table_input = ch_merged
        .map { m, table, fasta, lookup -> [ m, table, lookup ] }
        .join(TAXONOMY.out.taxonomy, by: 0)   // -> [ m, table, lookup, taxonomy ]

    ASV_TAXONOMY_TABLE(ch_table_input)
    ch_versions = ch_versions.mix(ASV_TAXONOMY_TABLE.out.versions)

    emit:
    asv_taxonomy_abundance = ASV_TAXONOMY_TABLE.out.abundance_table   // [ marker, {MARKER}.asv_taxonomy_abundance.tsv ]
    taxonomy_for_ecology   = ASV_TAXONOMY_TABLE.out.taxonomy_by_sequence // [ marker, {MARKER}.taxonomy_by_sequence.tsv ]
    versions               = ch_versions
}
