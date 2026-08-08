process HONEYPI_DADA2 {
    tag "honeypi_dada2"
    label 'process_high'

    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/pits/dada2", mode: 'copy'

    input:
    path(r1_files)   // collected R1 trimmed files
    path(r2_files)   // collected R2 trimmed files

    output:
    path 'ASVs.fasta',            emit: asvs
    path 'ASVs_dada2_counts.txt', emit: counts
    path 'dada2_stats.txt',       emit: stats
    path 'error_rates.pdf',       emit: error_plots, optional: true
    path 'versions.yml',          emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(dada2))

    trunc_r1  <- as.integer("${params.honeypi.dada2_trunc_len_r1}")
    trunc_r2  <- as.integer("${params.honeypi.dada2_trunc_len_r2}")
    ee_vals   <- as.numeric(strsplit("${params.honeypi.dada2_max_ee}", ",")[[1]])
    max_ee_r1 <- ee_vals[1]; max_ee_r2 <- ee_vals[2]
    n_cpus    <- as.integer("${task.cpus}")

    fnFs <- sort(list.files(".", pattern="val_1.fq.gz", full.names=TRUE))
    fnRs <- sort(list.files(".", pattern="val_2.fq.gz", full.names=TRUE))

    if (length(fnFs) == 0) stop("No trimmed R1 files found (pattern: val_1.fq.gz).")
    if (length(fnFs) != length(fnRs))
        stop(paste("R1/R2 count mismatch:", length(fnFs), "vs", length(fnRs)))

    sample_names <- sub("_R1_val_1[.]fq[.]gz\$", "", basename(fnFs))
    message("Processing ", length(sample_names), " samples: ",
            paste(head(sample_names, 5), collapse=", "),
            if (length(sample_names) > 5) "..." else "")

    # Filter and trim
    dir.create("filtered", showWarnings=FALSE)
    filtFs <- file.path("filtered", paste0(sample_names, "_F_filt.fq.gz"))
    filtRs <- file.path("filtered", paste0(sample_names, "_R_filt.fq.gz"))
    names(filtFs) <- sample_names
    names(filtRs) <- sample_names

    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
        truncLen  = c(trunc_r1, trunc_r2),
        maxN      = 0,
        maxEE     = c(max_ee_r1, max_ee_r2),
        truncQ    = 2,
        rm.phix   = TRUE,
        compress  = TRUE,
        multithread = n_cpus)

    write.table(data.frame(sample=rownames(out), out),
                "dada2_filter_stats.txt", sep="\\t", quote=FALSE, row.names=FALSE)
    message("Filter stats:\\n", paste(capture.output(print(out)), collapse="\\n"))

    keep <- out[, "reads.out"] >= 10
    if (!any(keep)) stop(
        "No samples retained after filtering. Check truncLen parameters. \\n",
        "Min reads required: 10. Stats:\\n",
        paste(capture.output(print(out)), collapse="\\n"))

    filtFs <- filtFs[keep]
    filtRs <- filtRs[keep]
    sample_names <- sample_names[keep]

    message("Retained ", sum(keep), " of ", length(keep), " samples after filtering.")

    # Error learning
    message("Learning error rates...")
    errF <- learnErrors(filtFs, multithread=n_cpus)
    errR <- learnErrors(filtRs, multithread=n_cpus)

    tryCatch({
        pdf("error_rates.pdf")
        plotErrors(errF, nominalQ=TRUE)
        dev.off()
    }, error=function(e) { message("Error plot failed: ", conditionMessage(e)) })

    # Dereplicate
    message("Dereplicating...")
    derepFs <- derepFastq(filtFs)
    derepRs <- derepFastq(filtRs)
    names(derepFs) <- sample_names
    names(derepRs) <- sample_names

    # DADA2 denoising
    message("Running DADA2 inference...")
    dadaFs <- dada(derepFs, err=errF, multithread=n_cpus)
    dadaRs <- dada(derepRs, err=errR, multithread=n_cpus)

    # Merge paired reads
    message("Merging paired reads...")
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                          maxMismatch=10, verbose=TRUE)

    # Sequence table
    seqtab <- makeSequenceTable(mergers)
    message("Amplicon length distribution:")
    print(table(nchar(getSequences(seqtab))))

    # Remove chimeras
    message("Removing chimeras...")
    seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                         multithread=n_cpus, verbose=TRUE)
    pct_removed <- round((1 - sum(seqtab_nochim) / sum(seqtab)) * 100, 1)
    message("Chimera removal: ", pct_removed, "% reads removed")

    if (ncol(seqtab_nochim) == 0)
        stop("No ASVs remaining after chimera removal.")

    # Assign ASV IDs
    seqs    <- colnames(seqtab_nochim)
    asv_ids <- sprintf("ASV_%010d", seq_along(seqs))

    # Write FASTA
    fasta_lines <- character(2L * length(seqs))
    for (i in seq_along(seqs)) {
        fasta_lines[2L * i - 1L] <- paste0(">", asv_ids[i])
        fasta_lines[2L * i]      <- seqs[i]
    }
    writeLines(fasta_lines, "ASVs.fasta")

    # Write count table (ASVs as rows, samples as columns)
    count_mat <- t(seqtab_nochim)
    rownames(count_mat) <- asv_ids
    write.table(data.frame(ASV_ID=rownames(count_mat), count_mat),
                "ASVs_dada2_counts.txt", sep="\\t", quote=FALSE, row.names=FALSE)

    # Summary stats
    n_input  <- sum(out[, "reads.in"])
    n_filter <- sum(out[keep, "reads.out"])
    n_final  <- sum(seqtab_nochim)
    stats_df <- data.frame(
        step  = c("input", "filtered", "non_chimeric"),
        reads = c(n_input, n_filter, n_final)
    )
    write.table(stats_df, "dada2_stats.txt", sep="\\t", quote=FALSE, row.names=FALSE)

    message("DADA2 complete: ", ncol(seqtab_nochim), " ASVs in ",
            nrow(seqtab_nochim), " samples")

    writeLines(c(
        '"HONEYPI_DADA2":',
        paste0('    dada2: ', packageVersion('dada2')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """

    stub:
    """
    printf '>ASV_0000000001\nACGTACGTACGTACGT\n>ASV_0000000002\nTGCATGCATGCATGCA\n' > ASVs.fasta
    printf 'ASV_ID\tsample-A\tsample-B\nASV_0000000001\t100\t80\nASV_0000000002\t50\t120\n' > ASVs_dada2_counts.txt
    printf 'step\treads\ninput\t10000\nfiltered\t9000\nnon_chimeric\t8500\n' > dada2_stats.txt
    printf '"HONEYPI_DADA2":\n    dada2: 1.28.0\n    R: 4.3.3\n' > versions.yml
    """
}
