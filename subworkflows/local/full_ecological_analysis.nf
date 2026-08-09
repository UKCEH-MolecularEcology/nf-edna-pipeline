include { ECOLOGY_ALPHA          } from '../../modules/local/ecology_alpha/main'
include { ECOLOGY_BETA           } from '../../modules/local/ecology_beta/main'
include { ECOLOGY_ORDINATION_FULL } from '../../modules/local/ecology_ordination_full/main'
include { ECOLOGY_DIFFERENTIAL   } from '../../modules/local/ecology_differential/main'
include { ECOLOGY_NETWORK        } from '../../modules/local/ecology_network/main'
include { ECOLOGY_INDICATORS     } from '../../modules/local/ecology_indicators/main'
include { ECOLOGY_ENVFIT         } from '../../modules/local/ecology_envfit/main'
include { ECOLOGY_MULTIMARKER    } from '../../modules/local/ecology_multimarker/main'
include { ECOLOGY_REPORT         } from '../../modules/local/ecology_report/main'

/*
 * FULL_ECOLOGICAL_ANALYSIS
 *
 * Comprehensive ecological analysis add-on running after the main eDNA pipeline.
 * Accepts the merged ASV tables and taxonomy tables from MERGE_ASV_TABLES and
 * TAXONOMY processes.
 *
 * Input:
 *   ch_ecology_input  — channel of [ marker, asv_table_tsv, taxonomy_tsv ]
 *   metadata          — path to sample metadata TSV (or assets/NO_FILE when absent)
 */

workflow FULL_ECOLOGICAL_ANALYSIS {

    take:
    ch_ecology_input   // [ marker, asv_table, taxonomy ]
    metadata           // path to metadata TSV, or assets/NO_FILE when absent

    main:
    ch_versions = Channel.empty()

    def meta_file = metadata
    def group_var = params.ecology_group_var ?: ''

    // ── Per-marker analyses (run in parallel for each marker) ────────────

    ECOLOGY_ALPHA(
        ch_ecology_input,
        meta_file,
        group_var
    )
    ch_versions = ch_versions.mix(ECOLOGY_ALPHA.out.versions.first())

    ECOLOGY_BETA(
        ch_ecology_input,
        meta_file,
        group_var
    )
    ch_versions = ch_versions.mix(ECOLOGY_BETA.out.versions.first())

    ECOLOGY_ORDINATION_FULL(
        ch_ecology_input,
        meta_file,
        group_var
    )
    ch_versions = ch_versions.mix(ECOLOGY_ORDINATION_FULL.out.versions.first())

    ECOLOGY_DIFFERENTIAL(
        ch_ecology_input,
        meta_file,
        group_var
    )
    ch_versions = ch_versions.mix(ECOLOGY_DIFFERENTIAL.out.versions.first())

    ECOLOGY_NETWORK(
        ch_ecology_input,
        params.ecology_min_prevalence ?: 0.3,
        params.ecology_cor_cutoff     ?: 0.6
    )
    ch_versions = ch_versions.mix(ECOLOGY_NETWORK.out.versions.first())

    ECOLOGY_INDICATORS(
        ch_ecology_input,
        meta_file,
        group_var
    )
    ch_versions = ch_versions.mix(ECOLOGY_INDICATORS.out.versions.first())

    ECOLOGY_ENVFIT(
        ch_ecology_input,
        meta_file,
        group_var
    )
    ch_versions = ch_versions.mix(ECOLOGY_ENVFIT.out.versions.first())

    // ── Cross-marker analysis (runs once all markers are done) ─────────
    ch_asv_tables_all  = ch_ecology_input.map { marker, asv, tax -> asv }.collect()
    ch_tax_tables_all  = ch_ecology_input.map { marker, asv, tax -> tax }.collect()
    ch_markers_list    = ch_ecology_input.map { marker, asv, tax -> marker }.collect()

    ECOLOGY_MULTIMARKER(
        ch_asv_tables_all,
        ch_tax_tables_all,
        ch_markers_list,
        meta_file
    )
    ch_versions = ch_versions.mix(ECOLOGY_MULTIMARKER.out.versions)

    // ── HTML report per marker ────────────────────────────────────────────
    // Join all per-marker result directories for the report
    ch_report_input = ECOLOGY_ALPHA.out.results
        .join(ECOLOGY_BETA.out.results,           by: 0)
        .join(ECOLOGY_ORDINATION_FULL.out.results, by: 0)
        .join(ECOLOGY_NETWORK.out.results,         by: 0)
        .join(ECOLOGY_INDICATORS.out.results,      by: 0)
        .join(ECOLOGY_ENVFIT.out.results,          by: 0)

    // Barrier: without this, a fast marker's report can render (and its
    // task can fail, e.g. an untriggered bug in the report script itself)
    // well before a slower marker has even finished DADA2, long before its
    // own ASV/taxonomy tables or ecology results exist. collect(flat: false)
    // waits for every marker's report inputs to arrive, keeping each one's
    // tuple intact (plain collect() would flatten all markers' items into
    // one scalar soup), then flatMap re-emits them individually -- so
    // reporting only starts once every marker is fully ready, not as each
    // one races in.
    ch_report_input
        .join(ch_ecology_input.map { m, a, t -> [m, a, t] }, by: 0)
        .collect(flat: false)
        .flatMap { it }
        .multiMap { m, alpha, beta, ord, net, ind, env, asv, tax ->
            dirs: [m, alpha, beta, ord, net, ind, env]
            asv:  asv
            tax:  tax
        }
        .set { rpt }

    ECOLOGY_REPORT(
        rpt.dirs,
        rpt.asv,
        rpt.tax,
        meta_file
    )
    ch_versions = ch_versions.mix(ECOLOGY_REPORT.out.versions.first())

    emit:
    alpha_results      = ECOLOGY_ALPHA.out.results
    beta_results       = ECOLOGY_BETA.out.results
    ordination_results = ECOLOGY_ORDINATION_FULL.out.results
    differential       = ECOLOGY_DIFFERENTIAL.out.results
    network_results    = ECOLOGY_NETWORK.out.results
    indicator_results  = ECOLOGY_INDICATORS.out.results
    envfit_results     = ECOLOGY_ENVFIT.out.results
    multimarker        = ECOLOGY_MULTIMARKER.out.results
    reports            = ECOLOGY_REPORT.out.report
    versions           = ch_versions
}
