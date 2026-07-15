process TAPIRS_BLAST_LCA {
    tag "${meta.id}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    input:
    tuple val(meta), path(blast_tsv)
    path taxdump_dir

    output:
    tuple val(meta), path('*.lca.tsv'), emit: lca
    path 'versions.yml',                emit: versions

    script:
    def prefix = "${meta.id}"
    def m = params.tapirs.mlca
    """
    taxdump_lineage.py \\
        --taxdump ${taxdump_dir} \\
        --mode blast \\
        --in ${blast_tsv} \\
        --out ${prefix}.blast.tax.tsv

    mlca.py \\
        --in ${prefix}.blast.tax.tsv \\
        --out ${prefix}.lca.tsv \\
        --identity ${m.identity} \\
        --coverage ${m.coverage} \\
        --bitscore ${m.bitscore} \\
        --majority ${m.majority} \\
        --min-hits ${m.min_hits}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.lca.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
