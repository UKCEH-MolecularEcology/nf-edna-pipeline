include { TAPIRS_BLAST             } from '../../modules/local/tapirs_blast/main'
include { TAPIRS_BLAST_LCA         } from '../../modules/local/tapirs_blast_lca/main'
include { LCA_TO_TAXONOMY_TABLE    } from '../../modules/local/lca_to_taxonomy_table/main'
include { TAPIRS_ASV_TAXONOMY_TABLE } from '../../modules/local/tapirs_asv_taxonomy_table/main'

// CO1 BLAST + majority-vote-LCA against nt (NCBI Nucleotide collection),
// reusing the same TAPIRS_BLAST / TAPIRS_BLAST_LCA modules 12S's Tapirs
// branch already uses -- pointed at the collective, post-merge CO1 ASV set
// (not a raw-read vsearch OTU-clustering front end, since BLAST classification
// only needs a query FASTA). Standalone output, mirrors 12S Tapirs: does NOT
// feed ecology. Uses the standard NCBI-taxid + taxdump lineage path (same as
// 12S), since nt hits carry real taxids -- unlike coidb, which needed the
// separate self-describing-header 'coidb' lineage mode.
workflow TAPIRS_BLAST_LCA_CO1 {

    take:
    ch_merged_co1   // [ marker, merged_table, merged_fasta, asv_lookup ] -- singleton, marker=='CO1'
    tapirs_opts     // params.tapirs_blast_lca_co1 { blast{db_dir,db_prefix,...}, mlca{...} }

    main:
    ch_versions = Channel.empty()

    ch_query = ch_merged_co1.map { marker, table, fasta, lookup -> [ [id: marker], fasta ] }

    TAPIRS_BLAST(
        ch_query,
        file(tapirs_opts.blast.db_dir),
        tapirs_opts.blast.db_prefix,
        tapirs_opts.blast
    )
    ch_versions = ch_versions.mix(TAPIRS_BLAST.out.versions)

    TAPIRS_BLAST_LCA(
        TAPIRS_BLAST.out.blast_tsv,
        file(tapirs_opts.taxdump),
        'blast',
        tapirs_opts.mlca
    )
    ch_versions = ch_versions.mix(TAPIRS_BLAST_LCA.out.versions)

    ch_lca_by_marker = TAPIRS_BLAST_LCA.out.lca.map { meta, tsv -> [ meta.id, tsv ] }

    LCA_TO_TAXONOMY_TABLE(ch_lca_by_marker)
    ch_versions = ch_versions.mix(LCA_TO_TAXONOMY_TABLE.out.versions)

    ch_table_input = ch_merged_co1
        .map { marker, table, fasta, lookup -> [ marker, table, lookup ] }
        .join(LCA_TO_TAXONOMY_TABLE.out.taxonomy, by: 0)   // -> [ marker, table, lookup, taxonomy ]

    TAPIRS_ASV_TAXONOMY_TABLE(ch_table_input)
    ch_versions = ch_versions.mix(TAPIRS_ASV_TAXONOMY_TABLE.out.versions)

    emit:
    asv_taxonomy_abundance = TAPIRS_ASV_TAXONOMY_TABLE.out.abundance_table
    taxonomy_by_sequence   = TAPIRS_ASV_TAXONOMY_TABLE.out.taxonomy_by_sequence
    versions               = ch_versions
}
