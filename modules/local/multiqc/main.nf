process MULTIQC {
    label 'process_single'

    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(multiqc_files, stageAs: 'multiqc_files/*')

    output:
    path '*multiqc_report.html', emit: report
    path '*_data',               emit: data
    path '*_plots',              emit: plots, optional: true
    path 'versions.yml',         emit: versions

    script:
    def args   = task.ext.args ?: ''
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    def title  = params.multiqc_title  ? "--title '${params.multiqc_title}'" : ''
    """
    multiqc \\
        --force \\
        ${config} \\
        ${title} \\
        ${args} \\
        multiqc_files/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    touch multiqc_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.21
    END_VERSIONS
    """
}
