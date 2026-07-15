process TAPIRS_VSEARCH_DENOISE_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0'

    input:
    tuple val(meta), path(derep)

    output:
    tuple val(meta), path('*.centroids.fasta'), emit: centroids
    path 'versions.yml',                        emit: versions

    script:
    def prefix = "${meta.id}"
    def v = params.tapirs.vsearch
    def method = params.tapirs.cluster_method
    if (method == 'cluster') {
        """
        vsearch --cluster_fast ${derep} \\
            --sizein --sizeout \\
            --query_cov ${v.query_cov} \\
            --id ${v.cluster_id} \\
            --strand both \\
            --centroids ${prefix}.centroids.fasta \\
            --fasta_width 0 \\
            --uc ${prefix}.cluster.tsv \\
            --threads ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//')
        END_VERSIONS
        """
    } else {
        """
        vsearch --cluster_unoise ${derep} \\
            --sizein --sizeout \\
            --minsize ${v.minsize} \\
            --unoise_alpha ${v.unoise_alpha} \\
            --id ${v.unoise_id} \\
            --centroids ${prefix}.centroids.fasta \\
            --fasta_width 0 \\
            --uc ${prefix}.denoise.tsv \\
            --threads ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.centroids.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: 2.27.0
    END_VERSIONS
    """
}
