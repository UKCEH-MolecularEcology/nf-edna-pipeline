process DADA2_FILTER {
    tag "${meta.id}_${marker}"
    label 'process_low'

    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/dada2/${marker}/filtered", mode: 'copy', pattern: '*.filter_stats.tsv'

    input:
    tuple val(run_id), val(meta), path(reads)
    val marker
    val trunc_len_f
    val trunc_len_r
    val max_ee_f
    val max_ee_r

    output:
    tuple val(run_id), val(meta), path('*.filtered.fastq.gz'), emit: filtered
    tuple val(run_id), val(meta), path('*.filter_stats.tsv'),  emit: stats
    path 'versions.yml',                                       emit: versions

    script:
    def prefix      = "${meta.id}_${marker}"
    def pe          = meta.single_end ? 'FALSE' : 'TRUE'
    def trunc_f_arg = trunc_len_f > 0 ? trunc_len_f : 0
    def trunc_r_arg = trunc_len_r > 0 ? trunc_len_r : 0
    """
    #!/usr/bin/env Rscript
    library(dada2)

    paired_end  <- as.logical("${pe}")
    trunc_len_f <- as.integer("${trunc_f_arg}")
    trunc_len_r <- as.integer("${trunc_r_arg}")
    max_ee_f    <- as.numeric("${max_ee_f}")
    max_ee_r    <- as.numeric("${max_ee_r}")

    all_files <- list.files(".", pattern = "\\\\.trimmed\\\\.fastq\\\\.gz\$", full.names = TRUE)

    # One sample per Nextflow task -- runs as a genuinely separate OS
    # process, so this is real parallelism across samples regardless of
    # multithread=TRUE/FALSE here. filterAndTrim()'s OWN multithread param
    # parallelises ACROSS FILES via mclapply (fork-based); with only one
    # file (pair) in this task there's nothing for it to fork over, so
    # multithread is left off entirely -- it would just add fork overhead
    # for zero benefit at this granularity.
    if (paired_end) {
        fwd_files <- sort(all_files[grepl("_R1", all_files)])
        rev_files <- sort(all_files[grepl("_R2", all_files)])
        fwd_filter <- gsub("trimmed", "filtered", fwd_files)
        rev_filter <- gsub("trimmed", "filtered", rev_files)

        out <- filterAndTrim(
            fwd_files, fwd_filter,
            rev_files, rev_filter,
            truncLen    = c(
                ifelse(trunc_len_f > 0, trunc_len_f, 0),
                ifelse(trunc_len_r > 0, trunc_len_r, 0)
            ),
            maxEE       = c(max_ee_f, max_ee_r),
            truncQ      = 2,
            minLen      = 50,
            rm.phix     = TRUE,
            compress    = TRUE
        )
    } else {
        fwd_filter <- gsub("trimmed", "filtered", all_files)
        out <- filterAndTrim(
            all_files, fwd_filter,
            maxEE       = max_ee_f,
            truncQ      = 2,
            minLen      = 50,
            rm.phix     = TRUE,
            compress    = TRUE
        )
    }

    write.table(data.frame(file = basename(rownames(out)), out, row.names = NULL),
                "${prefix}.filter_stats.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

    writeLines(
        c(
            paste0('"${task.process}":'),
            paste0('    dada2: ', packageVersion('dada2')),
            paste0('    R: ', R.version\$major, '.', R.version\$minor)
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = "${meta.id}_${marker}"
    """
    touch ${prefix}_R1.filtered.fastq.gz
    touch ${prefix}_R2.filtered.fastq.gz
    printf "file\treads.in\treads.out\n" > ${prefix}.filter_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dada2: 1.30.0
    END_VERSIONS
    """
}
