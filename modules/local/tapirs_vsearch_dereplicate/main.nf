process TAPIRS_VSEARCH_DEREPLICATE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.derep.fasta'), emit: derep
    path 'versions.yml',                    emit: versions

    script:
    def prefix = "${meta.id}"
    def v = params.tapirs.vsearch
    """
    vsearch --derep_fulllength ${fasta} \\
        --sizeout \\
        --minuniquesize ${v.minuniqsize} \\
        --output ${prefix}.derep.fasta \\
        --fasta_width 0 \\
        --uc ${prefix}.derep_clusters.tsv \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.derep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: 2.27.0
    END_VERSIONS
    """
}
