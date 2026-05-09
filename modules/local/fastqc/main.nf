process FASTQC {
    tag "${meta.id}_${meta.marker}"
    label 'process_medium'

    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir "${params.outdir}/fastqc/${meta.marker}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip'),  emit: zip
    path 'versions.yml',             emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.marker}"
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${args} \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}_${meta.marker}"
    """
    touch ${prefix}_fastqc.html
    touch ${prefix}_fastqc.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}
