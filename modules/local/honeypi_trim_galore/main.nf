process HONEYPI_TRIM_GALORE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

    publishDir "${params.outdir}/pits/trim_galore", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.fq.gz') ? "trimmed/${fn}" : "logs/${fn}" }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1_val_1.fq.gz"), path("${meta.id}_R2_val_2.fq.gz"), emit: reads
    path "${meta.id}_trim_stats.txt",                                                         emit: stats
    path 'versions.yml',                                                                       emit: versions

    script:
    """
    # Symlink to a fixed name so Trim Galore output is predictable,
    # regardless of what the input file is called after Nextflow staging.
    ln -sf "\$(readlink -f ${r1})" _r1_input.fastq.gz
    ln -sf "\$(readlink -f ${r2})" _r2_input.fastq.gz

    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --quality 20 \\
        --length 50 \\
        _r1_input.fastq.gz _r2_input.fastq.gz \\
        2>&1 | tee "${meta.id}_trim_stats.txt"

    # Trim Galore output: _r1_input_val_1.fq.gz  _r2_input_val_2.fq.gz
    mv _r1_input_val_1.fq.gz "${meta.id}_R1_val_1.fq.gz"
    mv _r2_input_val_2.fq.gz "${meta.id}_R2_val_2.fq.gz"

    trim_galore --version | grep 'Trim Galore' | sed 's/Trim Galore version //' | \\
        awk '{print "\\"HONEYPI_TRIM_GALORE\\":\\n    trim_galore: " \$1}' > versions.yml
    """

    stub:
    """
    touch "${meta.id}_R1_val_1.fq.gz" "${meta.id}_R2_val_2.fq.gz"
    touch "${meta.id}_trim_stats.txt"
    printf '"HONEYPI_TRIM_GALORE":\n    trim_galore: 0.6.10\n' > versions.yml
    """
}
