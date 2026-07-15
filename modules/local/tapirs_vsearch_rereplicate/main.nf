process TAPIRS_VSEARCH_REREPLICATE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0'

    publishDir "${params.outdir}/tapirs/09_rereplicated/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(nonchimeras)

    output:
    tuple val(meta), path('*.rerep.fasta'), emit: rerep
    path 'versions.yml',                    emit: versions

    script:
    def prefix = "${meta.id}"
    """
    vsearch --rereplicate ${nonchimeras} \\
        --fasta_width 0 \\
        --output ${prefix}.rerep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.rerep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: 2.27.0
    END_VERSIONS
    """
}
