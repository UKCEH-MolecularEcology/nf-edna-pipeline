process TAPIRS_FASTP {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    publishDir "${params.outdir}/tapirs/02_trimmed/${meta.id}", mode: 'copy', pattern: '*.fastp.json'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.forward.merged.fastq'), emit: fwd_merged
    tuple val(meta), path('*.trim.fastp.json'),      emit: trim_json
    path 'versions.yml',                             emit: versions

    script:
    def prefix = "${meta.id}"
    def f = params.tapirs.fastp
    """
    # Step 1: quality trim (no adapter trimming — Tapirs' own front-trim
    # substitutes for primer removal here)
    fastp \\
        --disable_adapter_trimming \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${prefix}.R1.trimmed.fastq \\
        --out2 ${prefix}.R2.trimmed.fastq \\
        --unpaired1 ${prefix}.R1.unpaired.fastq \\
        --unpaired2 ${prefix}.R2.unpaired.fastq \\
        --failed_out ${prefix}.trimmed.failed.fastq \\
        -j ${prefix}.trim.fastp.json \\
        -h ${prefix}.trim.fastp.html \\
        --qualified_quality_phred ${f.qual_phred} \\
        --unqualified_percent_limit ${f.unqualified_percent_limit} \\
        --average_qual ${f.qual_phred} \\
        --cut_tail \\
        --cut_window_size ${f.window_size} \\
        --cut_mean_quality ${f.qual_phred} \\
        --trim_poly_g \\
        --trim_poly_x \\
        --poly_g_min_len ${f.poly_g_min} \\
        --poly_x_min_len ${f.poly_x_min} \\
        --length_required ${f.len_required} \\
        --trim_front1 ${f.trim_front1} \\
        --trim_front2 ${f.trim_front2} \\
        --overlap_diff_percent_limit ${f.diff_percent_limit} \\
        --max_len1 ${f.max_len1} \\
        --max_len2 ${f.max_len2} \\
        --thread ${task.cpus}

    # Step 2: merge overlapping trimmed pairs
    fastp \\
        --disable_quality_filtering \\
        --disable_adapter_trimming \\
        --in1 ${prefix}.R1.trimmed.fastq \\
        --in2 ${prefix}.R2.trimmed.fastq \\
        --out1 ${prefix}.R1.unmerged.fastq \\
        --out2 ${prefix}.R2.unmerged.fastq \\
        --merge \\
        --merged_out ${prefix}.merged.fastq \\
        --overlap_len_require ${f.min_overlap} \\
        --overlap_diff_limit ${f.diff_limit} \\
        --overlap_diff_percent_limit ${f.diff_percent_limit} \\
        --length_limit ${f.length_limit} \\
        --length_required ${f.len_required} \\
        -j ${prefix}.merge.fastp.json \\
        -h ${prefix}.merge.fastp.html \\
        --correction \\
        --thread ${task.cpus}

    # Step 3: pool merged pairs + unpaired-but-passing R1 + unmerged R1 into
    # one forward-oriented file (mirrors merge_forward_reads in the Snakemake
    # pipeline, which just concatenates these 3 fastqs)
    cat ${prefix}.merged.fastq ${prefix}.R1.unpaired.fastq ${prefix}.R1.unmerged.fastq \\
        > ${prefix}.forward.merged.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.forward.merged.fastq
    touch ${prefix}.trim.fastp.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 0.23.4
    END_VERSIONS
    """
}
