process TAPIRS_VSEARCH_CHIMERA {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0'

    input:
    tuple val(meta), path(centroids)
    path dechim_db   // only read when params.tapirs.chimera_detection == 'ref'

    output:
    tuple val(meta), path('*.nc.fasta'),      emit: nonchimeras
    tuple val(meta), path('*.chimera.fasta'), emit: chimeras
    path 'versions.yml',                      emit: versions

    script:
    def prefix = "${meta.id}"
    def v = params.tapirs.vsearch
    if (params.tapirs.chimera_detection == 'ref') {
        """
        vsearch --uchime_ref ${centroids} \\
            --db ${dechim_db} \\
            --chimeras ${prefix}.chimera.fasta \\
            --borderline ${prefix}.chimera.fasta \\
            --mindiffs ${v.mindiffs} \\
            --mindiv ${v.mindiv} \\
            --fasta_width 0 \\
            --nonchimeras ${prefix}.nc.fasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//')
        END_VERSIONS
        """
    } else {
        """
        vsearch --uchime3_denovo ${centroids} \\
            --abskew ${v.abskew} \\
            --chimeras ${prefix}.chimera.fasta \\
            --borderline ${prefix}.chimera.fasta \\
            --fasta_width 0 \\
            --nonchimeras ${prefix}.nc.fasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.nc.fasta
    touch ${prefix}.chimera.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: 2.27.0
    END_VERSIONS
    """
}
