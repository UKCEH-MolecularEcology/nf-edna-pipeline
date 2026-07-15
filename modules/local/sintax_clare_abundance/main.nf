process SINTAX_CLARE_ABUNDANCE {
    tag "${marker}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    input:
    tuple val(marker), path(clare_abundance_csv)

    output:
    tuple val(marker), path('*.ncl_matrix_raw.csv'), emit: ncl_matrix
    path 'versions.yml',                             emit: versions

    script:
    def prefix = "${marker}"
    """
    sintax_clare_abundance.py \\
        --in ${clare_abundance_csv} \\
        --out ${prefix}.ncl_matrix_raw.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${marker}"
    """
    touch ${prefix}.ncl_matrix_raw.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
