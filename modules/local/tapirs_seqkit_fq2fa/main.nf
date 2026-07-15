process TAPIRS_SEQKIT_FQ2FA {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path 'versions.yml',              emit: versions

    script:
    def prefix = "${meta.id}"
    """
    seqkit fq2fa ${fastq} -o ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.8.2
    END_VERSIONS
    """
}
