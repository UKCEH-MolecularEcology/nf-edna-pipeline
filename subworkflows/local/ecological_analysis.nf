include { ECOLOGY_DIVERSITY  } from '../../modules/local/ecology_diversity/main'
include { ECOLOGY_ORDINATION } from '../../modules/local/ecology_ordination/main'
include { ECOLOGY_BARPLOT    } from '../../modules/local/ecology_barplot/main'

workflow ECOLOGICAL_ANALYSIS {

    take:
    ch_input    // [ marker, asv_table, taxonomy_table ]
    metadata    // file or [] if none

    main:

    // Alpha & beta diversity, statistical tests
    ECOLOGY_DIVERSITY(ch_input, metadata)

    // Ordination: PCoA, NMDS
    ECOLOGY_ORDINATION(ch_input, metadata)

    // Taxonomic composition barplots
    ECOLOGY_BARPLOT(ch_input, metadata)

    emit:
    diversity_results  = ECOLOGY_DIVERSITY.out.results
    ordination_results = ECOLOGY_ORDINATION.out.results
    barplot_results    = ECOLOGY_BARPLOT.out.results
}
