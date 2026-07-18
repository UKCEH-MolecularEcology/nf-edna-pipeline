include { TAPIRS_FASTP                  } from '../../modules/local/tapirs_fastp/main'
include { TAPIRS_SEQKIT_FQ2FA           } from '../../modules/local/tapirs_seqkit_fq2fa/main'
include { TAPIRS_VSEARCH_DEREPLICATE    } from '../../modules/local/tapirs_vsearch_dereplicate/main'
include { TAPIRS_VSEARCH_DENOISE_CLUSTER} from '../../modules/local/tapirs_vsearch_denoise_cluster/main'
include { TAPIRS_VSEARCH_CHIMERA        } from '../../modules/local/tapirs_vsearch_chimera/main'
include { TAPIRS_VSEARCH_REREPLICATE    } from '../../modules/local/tapirs_vsearch_rereplicate/main'
include { TAPIRS_BLAST                  } from '../../modules/local/tapirs_blast/main'
include { TAPIRS_BLAST_LCA              } from '../../modules/local/tapirs_blast_lca/main'
include { TAPIRS_KRAKEN2                } from '../../modules/local/tapirs_kraken2/main'
include { TAPIRS_KRAKEN2_TAXONOMY       } from '../../modules/local/tapirs_kraken2_taxonomy/main'
include { TAPIRS_OTU_TABLE as TAPIRS_OTU_TABLE_BLAST   } from '../../modules/local/tapirs_otu_table/main'
include { TAPIRS_OTU_TABLE as TAPIRS_OTU_TABLE_KRAKEN2 } from '../../modules/local/tapirs_otu_table/main'

// The standalone 12S "Tapirs" branch (fastp/vsearch/BLAST/Kraken2/MLCA),
// ported from 12S-edna-dada2-tapirs-workflow. Fed RAW reads (not cutadapt
// output) since fastp's own front-trim substitutes for primer removal here;
// runs independently of the DADA2/SINTAX branch and does not feed ecology.
workflow TAPIRS_12S {

    take:
    ch_reads   // [ meta, reads ] — 12S samples only, untrimmed

    main:
    ch_versions = Channel.empty()
    t = params.tapirs

    def missing = []
    if (!t.blast.db && t.analysis_method in ['blast', 'both'])       missing << 'tapirs.blast.db'
    if (!t.kraken2.db && t.analysis_method in ['kraken2', 'both'])   missing << 'tapirs.kraken2.db'
    if (!t.taxdump)                                                   missing << 'tapirs.taxdump'
    if (t.chimera_detection == 'ref' && !t.dechim_blast_db)           missing << 'tapirs.dechim_blast_db'
    if (!t.experiment)                                                 missing << 'tapirs.experiment'
    if (missing) {
        error "params.tapirs.enabled is true but required path(s) are not set: ${missing.join(', ')}"
    }

    ch_dechim_db = Channel.value(
        t.chimera_detection == 'ref' ? file(t.dechim_blast_db) : file("${projectDir}/assets/NO_FILE")
    )
    ch_taxdump_dir = Channel.value(file(t.taxdump))

    // 1. fastp trim + merge → forward-oriented pooled FASTQ
    TAPIRS_FASTP(ch_reads)
    ch_versions = ch_versions.mix(TAPIRS_FASTP.out.versions.first())

    // 2. FASTQ → FASTA
    TAPIRS_SEQKIT_FQ2FA(TAPIRS_FASTP.out.fwd_merged)
    ch_versions = ch_versions.mix(TAPIRS_SEQKIT_FQ2FA.out.versions.first())

    // 3. Dereplicate
    TAPIRS_VSEARCH_DEREPLICATE(TAPIRS_SEQKIT_FQ2FA.out.fasta)
    ch_versions = ch_versions.mix(TAPIRS_VSEARCH_DEREPLICATE.out.versions.first())

    // 4. Denoise (UNOISE3) or cluster
    TAPIRS_VSEARCH_DENOISE_CLUSTER(TAPIRS_VSEARCH_DEREPLICATE.out.derep)
    ch_versions = ch_versions.mix(TAPIRS_VSEARCH_DENOISE_CLUSTER.out.versions.first())

    // 5. Chimera removal (denovo or ref)
    TAPIRS_VSEARCH_CHIMERA(TAPIRS_VSEARCH_DENOISE_CLUSTER.out.centroids, ch_dechim_db)
    ch_versions = ch_versions.mix(TAPIRS_VSEARCH_CHIMERA.out.versions.first())

    // 6. Rereplicate — always runs: its output is a required read-count
    // denominator for the BLAST-flavored OTU table too, not just Kraken2's.
    TAPIRS_VSEARCH_REREPLICATE(TAPIRS_VSEARCH_CHIMERA.out.nonchimeras)
    ch_versions = ch_versions.mix(TAPIRS_VSEARCH_REREPLICATE.out.versions.first())

    ch_rerep_all = TAPIRS_VSEARCH_REREPLICATE.out.rerep.map { meta, fa -> fa }.collect()

    ch_otu_blast    = Channel.empty()
    ch_otu_kraken2  = Channel.empty()

    // BLAST → taxdump lineage → majority-vote LCA → OTU table
    if (t.analysis_method in ['blast', 'both']) {
        def blast_db_path   = file(t.blast.db)
        ch_blast_db_dir    = Channel.value(blast_db_path.getParent())
        ch_blast_db_prefix = Channel.value(blast_db_path.getName())

        TAPIRS_BLAST(TAPIRS_VSEARCH_CHIMERA.out.nonchimeras, ch_blast_db_dir, ch_blast_db_prefix, t.blast)
        ch_versions = ch_versions.mix(TAPIRS_BLAST.out.versions.first())

        TAPIRS_BLAST_LCA(TAPIRS_BLAST.out.blast_tsv, ch_taxdump_dir, 'blast', t.mlca)
        ch_versions = ch_versions.mix(TAPIRS_BLAST_LCA.out.versions.first())

        ch_lca_all = TAPIRS_BLAST_LCA.out.lca.map { meta, tsv -> tsv }.collect()

        TAPIRS_OTU_TABLE_BLAST(
            ch_lca_all,
            ch_rerep_all,
            "${t.experiment}_blast${t.mlca.identity}_${t.cluster_method}",
            t.lowest_taxonomic_rank,
            t.highest_taxonomic_rank
        )
        ch_versions = ch_versions.mix(TAPIRS_OTU_TABLE_BLAST.out.versions)
        ch_otu_blast = TAPIRS_OTU_TABLE_BLAST.out.otu_table
    }

    // Kraken2 → taxdump lineage → OTU table
    if (t.analysis_method in ['kraken2', 'both']) {
        ch_kraken2_db_dir = Channel.value(file(t.kraken2.db))

        TAPIRS_KRAKEN2(TAPIRS_VSEARCH_REREPLICATE.out.rerep, ch_kraken2_db_dir)
        ch_versions = ch_versions.mix(TAPIRS_KRAKEN2.out.versions.first())

        TAPIRS_KRAKEN2_TAXONOMY(TAPIRS_KRAKEN2.out.output, ch_taxdump_dir)
        ch_versions = ch_versions.mix(TAPIRS_KRAKEN2_TAXONOMY.out.versions.first())

        ch_krk_tax_all = TAPIRS_KRAKEN2_TAXONOMY.out.tax.map { meta, tsv -> tsv }.collect()

        def conf_decimal = t.kraken2.confidence.toString().split('\\.')[1]
        TAPIRS_OTU_TABLE_KRAKEN2(
            ch_krk_tax_all,
            ch_rerep_all,
            "${t.experiment}_kraken2_conf${conf_decimal}_${t.cluster_method}",
            t.lowest_taxonomic_rank,
            t.highest_taxonomic_rank
        )
        ch_versions = ch_versions.mix(TAPIRS_OTU_TABLE_KRAKEN2.out.versions)
        ch_otu_kraken2 = TAPIRS_OTU_TABLE_KRAKEN2.out.otu_table
    }

    emit:
    otu_table_blast   = ch_otu_blast
    otu_table_kraken2 = ch_otu_kraken2
    versions          = ch_versions
}
