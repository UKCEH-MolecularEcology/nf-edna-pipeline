process HONEYPI_SUMMARY {
    tag "summary"
    label 'process_low'

    publishDir "${params.outdir}/pits/honeypi_output", mode: 'copy'

    input:
    path asvs_fasta,       stageAs: 'in_asvs.fasta'
    path counts_merged,    stageAs: 'in_counts_merged.txt'
    path counts_filtered,  stageAs: 'in_counts_filtered.txt'
    path taxonomy,         stageAs: 'in_taxonomy.txt'
    path dada2_stats,      stageAs: 'in_dada2_stats.txt'
    path error_rates,      stageAs: 'in_error_rates.pdf'

    output:
    path 'ASVs.fasta'
    path 'ASVs_counts_merged.txt'
    path 'ASVs_counts_filtered.txt'
    path 'ASVs_taxonomy.txt'
    path 'dada2_stats.txt'
    path 'error_rates.pdf', optional: true

    script:
    """
    cp in_asvs.fasta          ASVs.fasta
    cp in_counts_merged.txt   ASVs_counts_merged.txt
    cp in_counts_filtered.txt ASVs_counts_filtered.txt
    cp in_taxonomy.txt        ASVs_taxonomy.txt
    cp in_dada2_stats.txt     dada2_stats.txt
    [ "${error_rates}" != "NO_FILE" ] && cp in_error_rates.pdf error_rates.pdf || true
    """

    stub:
    """
    touch ASVs.fasta ASVs_counts_merged.txt ASVs_counts_filtered.txt
    touch ASVs_taxonomy.txt dada2_stats.txt
    """
}
