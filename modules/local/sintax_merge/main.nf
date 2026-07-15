process SINTAX_MERGE {
    tag "${marker}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/sintax/${marker}", mode: 'copy'

    input:
    tuple val(marker), path(asv_lookup), val(db_names), path(parsed_tsvs)

    output:
    tuple val(marker), path('*.asv_taxonomy_compare.tsv'),    emit: compare
    tuple val(marker), path('*.sintax_database_summary.tsv'), emit: summary
    path 'versions.yml',                                       emit: versions

    script:
    def prefix = "${marker}"
    def pairs = [db_names, parsed_tsvs].transpose().collect { db, f -> "${db}=${f}" }.join(' ')
    """
    sintax_merge.py \\
        --asv-lookup ${asv_lookup} \\
        --parsed ${pairs} \\
        --out-compare ${prefix}.asv_taxonomy_compare.tsv \\
        --out-summary ${prefix}.sintax_database_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${marker}"
    """
    touch ${prefix}.asv_taxonomy_compare.tsv
    touch ${prefix}.sintax_database_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
