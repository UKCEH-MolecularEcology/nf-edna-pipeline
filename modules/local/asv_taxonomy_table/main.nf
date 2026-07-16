process ASV_TAXONOMY_TABLE {
    tag "${marker}"
    label 'process_single'

    // Generic python3+pandas image (built for the 12S Tapirs/SINTAX branch,
    // but has no 12S-specific content) -- reused here to avoid a new
    // container just for a lightweight TSV join.
    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/asv_taxonomy/${marker}", mode: 'copy'

    input:
    tuple val(marker), path(merged_table), path(asv_lookup), path(taxonomy)

    output:
    tuple val(marker), path('*.asv_taxonomy_abundance.tsv'), emit: abundance_table
    tuple val(marker), path('*.taxonomy_by_sequence.tsv'),   emit: taxonomy_by_sequence
    path 'versions.yml',                                     emit: versions

    script:
    def prefix = "${marker}"
    """
    build_taxonomy_table.py \\
        --asv-lookup ${asv_lookup} \\
        --merged-table ${merged_table} \\
        --taxonomy ${taxonomy} \\
        --out-abundance ${prefix}.asv_taxonomy_abundance.tsv \\
        --out-by-sequence ${prefix}.taxonomy_by_sequence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${marker}"
    """
    touch ${prefix}.asv_taxonomy_abundance.tsv
    touch ${prefix}.taxonomy_by_sequence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
