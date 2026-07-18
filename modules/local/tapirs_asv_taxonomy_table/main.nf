process TAPIRS_ASV_TAXONOMY_TABLE {
    tag "${marker}"
    label 'process_single'

    // Same join script/container as the generic ASV_TAXONOMY_TABLE module --
    // build_taxonomy_table.py only cares about column shape, not where the
    // taxonomy came from. A separate process (rather than reusing
    // ASV_TAXONOMY_TABLE directly) exists purely so this can publish
    // alongside the other Tapirs BLAST+LCA output instead of colliding with
    // the DADA2-based asv_taxonomy/{MARKER}/ files for the same marker.
    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/tapirs_blast/${marker}", mode: 'copy'

    input:
    tuple val(marker), path(merged_table), path(asv_lookup), path(taxonomy)

    output:
    tuple val(marker), path('*.asv_taxonomy_abundance.tsv'), emit: abundance_table
    tuple val(marker), path('*.taxonomy_by_sequence.tsv'),   emit: taxonomy_by_sequence
    path 'versions.yml',                                     emit: versions

    script:
    def prefix = "${marker}_blast"
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
    def prefix = "${marker}_blast"
    """
    touch ${prefix}.asv_taxonomy_abundance.tsv
    touch ${prefix}.taxonomy_by_sequence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
