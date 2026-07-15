process TAPIRS_KRAKEN2_TAXONOMY {
    tag "${meta.id}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    input:
    tuple val(meta), path(krk_output)
    path taxdump_dir

    output:
    tuple val(meta), path('*.krk.tax.tsv'), emit: tax
    path 'versions.yml',                    emit: versions

    script:
    def prefix = "${meta.id}"
    """
    taxdump_lineage.py \\
        --taxdump ${taxdump_dir} \\
        --mode kraken2 \\
        --in ${krk_output} \\
        --out ${prefix}.krk.tax.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.krk.tax.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
