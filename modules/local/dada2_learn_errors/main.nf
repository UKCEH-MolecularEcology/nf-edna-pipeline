process DADA2_LEARN_ERRORS {
    tag "${run_id}_${marker}"
    label 'process_high'

    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/dada2/${marker}/error_models", mode: 'copy'

    input:
    tuple val(run_id), val(metas), path(filtered_reads)
    val marker

    output:
    tuple val(run_id), path('error_model_fwd.rds'), path('error_model_rev.rds'), emit: error_model
    path '*.png',                                                                emit: plots
    path 'versions.yml',                                                        emit: versions

    script:
    def pe = metas[0].single_end ? 'FALSE' : 'TRUE'
    """
    #!/usr/bin/env Rscript
    library(dada2)

    paired_end <- as.logical("${pe}")

    # Samples were already filtered per-sample by DADA2_FILTER (each its own
    # Nextflow task -- real parallelism across samples, unlike relying on
    # filterAndTrim's own multithread param, which forks via mclapply and
    # has nothing to fork over once filtering happens one sample at a time).
    all_filtered <- list.files(".", pattern = "\\\\.filtered\\\\.fastq\\\\.gz\$", full.names = TRUE)

    if (paired_end) {
        fwd_filter <- sort(all_filtered[grepl("_R1", all_filtered)])
        rev_filter <- sort(all_filtered[grepl("_R2", all_filtered)])

        fwd_passed <- fwd_filter[file.exists(fwd_filter) & file.info(fwd_filter)\$size > 0]
        rev_passed <- rev_filter[file.exists(rev_filter) & file.info(rev_filter)\$size > 0]

        # multithread=TRUE: matches the reference DADA2 big-data workflow's
        # own learnErrors() call exactly (no randomize -- that isn't in the
        # reference call either, and defaults to FALSE).
        err_fwd <- learnErrors(fwd_passed, multithread = ${task.cpus})
        err_rev <- learnErrors(rev_passed, multithread = ${task.cpus})

    } else {
        fwd_passed <- all_filtered[file.exists(all_filtered) & file.info(all_filtered)\$size > 0]
        err_fwd <- learnErrors(fwd_passed, multithread = ${task.cpus})
        err_rev <- err_fwd  # Placeholder for SE
    }

    saveRDS(err_fwd, "error_model_fwd.rds")
    saveRDS(err_rev, "error_model_rev.rds")

    # Plot error models
    png("error_model_fwd.png", width = 1200, height = 900)
    plotErrors(err_fwd, nominalQ = TRUE)
    dev.off()

    if (paired_end) {
        png("error_model_rev.png", width = 1200, height = 900)
        plotErrors(err_rev, nominalQ = TRUE)
        dev.off()
    }

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
    """
    touch error_model_fwd.rds
    touch error_model_rev.rds
    touch error_models.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dada2: 1.30.0
    END_VERSIONS
    """

}
