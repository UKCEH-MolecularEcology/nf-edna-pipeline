process FASTQC {
    tag "${meta.id}_${meta.marker}"
    label 'process_medium'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

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

    # FastQC names its outputs after the input file basename, which doesn't
    # include marker -- when the same sample (same DNA extract) is run
    # through multiple markers, MultiQC's collected file list ends up with
    # identically-named zips from different markers and errors out with an
    # "input file name collision". Prefix with the marker to disambiguate.
    for f in *_fastqc.html *_fastqc.zip; do
        [ -e "\$f" ] || continue
        mv "\$f" "${meta.marker}_\$f"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}_${meta.marker}"
    """
    touch ${meta.marker}_${prefix}_fastqc.html
    touch ${meta.marker}_${prefix}_fastqc.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}
