process DADA2_LEARN_ERRORS {
    tag "${run_id}_${marker}"
    label 'process_high'

    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/dada2/${marker}/error_models", mode: 'copy'

    input:
    tuple val(run_id), val(metas), path(reads)
    val marker
    val trunc_len_f
    val trunc_len_r
    val max_ee_f
    val max_ee_r

    output:
    tuple val(run_id), path('error_model_fwd.rds'), path('error_model_rev.rds'), emit: error_model
    path '*.png',                                                                 emit: plots
    path 'versions.yml',                                                         emit: versions

    script:
    def pe          = metas[0].single_end ? 'FALSE' : 'TRUE'
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

    # Collect all trimmed read files for this run
    all_files <- list.files(".", pattern = "\\\\.trimmed\\\\.fastq\\\\.gz\$", full.names = TRUE)

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
            compress    = TRUE,
            multithread = ${task.cpus}
        )
        message("Filter results (fwd):"); print(out)

        # Keep only samples with reads after filtering
        fwd_passed <- fwd_filter[file.exists(fwd_filter) & file.info(fwd_filter)\$size > 0]
        rev_passed <- rev_filter[file.exists(rev_filter) & file.info(rev_filter)\$size > 0]

        err_fwd <- learnErrors(fwd_passed, multithread = ${task.cpus}, randomize = TRUE)
        err_rev <- learnErrors(rev_passed, multithread = ${task.cpus}, randomize = TRUE)

    } else {
        fwd_filter <- gsub("trimmed", "filtered", all_files)
        out <- filterAndTrim(
            all_files, fwd_filter,
            maxEE       = max_ee_f,
            truncQ      = 2,
            minLen      = 50,
            rm.phix     = TRUE,
            compress    = TRUE,
            multithread = ${task.cpus}
        )
        fwd_passed <- fwd_filter[file.exists(fwd_filter) & file.info(fwd_filter)\$size > 0]
        err_fwd <- learnErrors(fwd_passed, multithread = ${task.cpus}, randomize = TRUE)
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
