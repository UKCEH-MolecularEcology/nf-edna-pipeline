process LCA_TO_TAXONOMY_TABLE {
    tag "${marker}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    input:
    tuple val(marker), path(lca_tsv)

    output:
    tuple val(marker), path('*.taxonomy.tsv'), emit: taxonomy
    path 'versions.yml',                       emit: versions

    script:
    def prefix = "${marker}"
    """
    lca_to_taxonomy_table.py \\
        --in ${lca_tsv} \\
        --out ${prefix}.taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${marker}"
    """
    touch ${prefix}.taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
