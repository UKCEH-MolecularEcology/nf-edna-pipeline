include { TAXONOMY               } from '../../modules/local/taxonomy/main'
include { MERGE_TAXONOMY_CHUNKS  } from '../../modules/local/merge_taxonomy_chunks/main'
include { ASV_TAXONOMY_TABLE     } from '../../modules/local/asv_taxonomy_table/main'

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

    ch_fasta = ch_merged.map { m, table, fasta, lookup -> [ m, fasta ] }

    // assignTaxonomy() (dada2 method)'s own multithread param forks via
    // mclapply -- the same mechanism already confirmed to give zero
    // speedup in this Singularity container for filterAndTrim (see
    // DADA2_FILTER). Split the collective ASV set into small chunks and
    // classify each as its own Nextflow task instead, so real parallelism
    // comes from many concurrent low-cpu tasks rather than one large task
    // stuck on an effectively single core. rdp/blast aren't affected
    // (rdp_classifier/blastn's own threading genuinely works) so they stay
    // as one task on the whole merged FASTA.
    ch_tax_input = marker_params.tax_method == 'dada2'
        ? ch_fasta.splitFasta(by: params.tax_chunk_size, file: true, elem: 1)
        : ch_fasta

    TAXONOMY(
        ch_tax_input,
        marker_params.tax_db,
        marker_params.tax_db_type,
        marker_params.tax_method,
        marker_params.addspecies_db
    )
    ch_versions = ch_versions.mix(TAXONOMY.out.versions)

    if (marker_params.tax_method == 'dada2') {
        MERGE_TAXONOMY_CHUNKS(marker, TAXONOMY.out.taxonomy.map { m, tsv -> tsv }.collect())
        ch_taxonomy = MERGE_TAXONOMY_CHUNKS.out.taxonomy
    } else {
        ch_taxonomy = TAXONOMY.out.taxonomy
    }

    ch_table_input = ch_merged
        .map { m, table, fasta, lookup -> [ m, table, lookup ] }
        .join(ch_taxonomy, by: 0)   // -> [ m, table, lookup, taxonomy ]

    ASV_TAXONOMY_TABLE(ch_table_input)
    ch_versions = ch_versions.mix(ASV_TAXONOMY_TABLE.out.versions)

    emit:
    asv_taxonomy_abundance = ASV_TAXONOMY_TABLE.out.abundance_table   // [ marker, {MARKER}.asv_taxonomy_abundance.tsv ]
    taxonomy_for_ecology   = ASV_TAXONOMY_TABLE.out.taxonomy_by_sequence // [ marker, {MARKER}.taxonomy_by_sequence.tsv ]
    versions               = ch_versions
}
