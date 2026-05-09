process VSEARCH_CHIMERA {
    tag "${meta.id}_${marker}"
    label 'process_medium'

    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0'

    publishDir "${params.outdir}/chimera_check/${marker}", mode: 'copy'

    input:
    tuple val(meta), path(fasta)
    val marker

    output:
    tuple val(meta), path('*.nonchimeras.fasta'), emit: nonchimeras
    tuple val(meta), path('*.chimera_ids.txt'),   emit: nonchimera_ids
    tuple val(meta), path('*.chimera_stats.txt'), emit: stats
    path 'versions.yml',                           emit: versions

    script:
    def prefix = "${meta.id}_${marker}"
    def id     = params.vsearch_id ?: 1.0
    """
    # De novo chimera detection (uchime3 algorithm)
    vsearch \\
        --uchime3_denovo ${fasta} \\
        --sizein \\
        --sizeout \\
        --fasta_width 0 \\
        --qmask none \\
        --nonchimeras ${prefix}.nonchimeras.fasta \\
        --chimeras ${prefix}.chimeras.fasta \\
        --uchimealns ${prefix}.chimera_alignments.txt \\
        --uchimeout ${prefix}.chimera_stats.txt \\
        --threads ${task.cpus} \\
        2>&1 | tee ${prefix}.vsearch.log

    # Extract non-chimeric ASV IDs for downstream table filtering
    grep "^>" ${prefix}.nonchimeras.fasta | sed 's/^>//' | cut -d';' -f1 \\
        > ${prefix}.chimera_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -1 | sed 's/vsearch v//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}_${meta.marker}"
    """
    touch ${prefix}.nonchimeras.fasta
    touch ${prefix}.chimera_ids.txt
    touch ${prefix}.chimera_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: 2.26.0
    END_VERSIONS
    """
}
