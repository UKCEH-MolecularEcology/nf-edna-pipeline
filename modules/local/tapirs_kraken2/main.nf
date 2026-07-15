process TAPIRS_KRAKEN2 {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0'

    input:
    tuple val(meta), path(rerep_fasta)
    path kraken2_db_dir

    output:
    tuple val(meta), path('*.krk'), emit: output
    tuple val(meta), path('*.txt'), emit: report
    path 'versions.yml',            emit: versions

    script:
    def prefix = "${meta.id}"
    def k = params.tapirs.kraken2
    """
    if [ -s ${rerep_fasta} ]; then
        kraken2 --db ${kraken2_db_dir} ${rerep_fasta} \\
            --threads ${task.cpus} \\
            --confidence ${k.confidence} \\
            --report ${prefix}.txt \\
            --output ${prefix}.krk
    else
        touch ${prefix}.txt
        touch ${prefix}.krk
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -1 | sed 's/^Kraken version //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.krk
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: 2.1.3
    END_VERSIONS
    """
}
