process DADA2_DENOISE {
    tag "${meta.id}_${marker}"
    label 'process_high'

    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/dada2/${marker}/asv_tables", mode: 'copy'

    input:
    tuple val(run_id), val(meta), path(reads), path(err_fwd), path(err_rev)
    val marker
    val trunc_len_f
    val trunc_len_r
    val max_ee_f
    val max_ee_r
    val pool
    val min_length
    val max_length

    output:
    tuple val(meta), path('*.asv_table.rds'),   emit: asv_table
    tuple val(meta), path('*.asv_seqs.fasta'),  emit: asv_seqs
    tuple val(meta), path('*.read_stats.tsv'),  emit: read_stats
    path 'versions.yml',                         emit: versions

    script:
    def prefix      = "${meta.id}_${marker}"
    def pe          = meta.single_end ? 'FALSE' : 'TRUE'
    def pool_arg    = pool == true ? 'TRUE' : pool == 'pseudo' ? '"pseudo"' : 'FALSE'
    def trunc_f_arg = trunc_len_f > 0 ? trunc_len_f : 0
    def trunc_r_arg = trunc_len_r > 0 ? trunc_len_r : 0
    """
    #!/usr/bin/env Rscript
    library(dada2)

    sample_id   <- "${meta.id}"
    marker      <- "${marker}"
    paired_end  <- as.logical("${pe}")
    trunc_len_f <- as.integer("${trunc_f_arg}")
    trunc_len_r <- as.integer("${trunc_r_arg}")
    max_ee_f    <- as.numeric("${max_ee_f}")
    max_ee_r    <- as.numeric("${max_ee_r}")
    pool        <- ${pool_arg}
    min_length  <- as.integer("${min_length}")
    max_length  <- as.integer("${max_length}")

    err_fwd <- readRDS("${err_fwd}")
    err_rev <- readRDS("${err_rev}")

    all_files <- list.files(".", pattern = "\\\\.trimmed\\\\.fastq\\\\.gz\$", full.names = TRUE)

    if (paired_end) {
        fwd_files <- sort(all_files[grepl("_R1", all_files)])
        rev_files <- sort(all_files[grepl("_R2", all_files)])

        fwd_filter <- gsub("trimmed", "filtered", fwd_files)
        rev_filter <- gsub("trimmed", "filtered", rev_files)

        # Per-read minLen: use truncation length when truncating (reads will be
        # exactly truncLen after filtering). min_length/max_length describe the
        # final MERGED amplicon length (applied below, after merging) -- they
        # are the wrong floor for an individual pre-merge R1/R2 read whenever
        # the amplicon is longer than one read (the untruncated/variable-length
        # case, trunc_len == 0), so fall back to DADA2's own conventional
        # default (20) instead of min_length there.
        min_per_read <- min(
            ifelse(trunc_len_f > 0, trunc_len_f, 20L),
            ifelse(trunc_len_r > 0, trunc_len_r, 20L)
        )

        # Re-filter for this sample
        out_filter <- filterAndTrim(
            fwd_files, fwd_filter,
            rev_files, rev_filter,
            truncLen    = c(
                ifelse(trunc_len_f > 0, trunc_len_f, 0),
                ifelse(trunc_len_r > 0, trunc_len_r, 0)
            ),
            maxEE       = c(max_ee_f, max_ee_r),
            truncQ      = 2,
            minLen      = min_per_read,
            rm.phix     = TRUE,
            compress    = TRUE,
            multithread = ${task.cpus}
        )

        fwd_ok <- fwd_filter[file.exists(fwd_filter) & file.info(fwd_filter)\$size > 0]
        rev_ok <- rev_filter[file.exists(rev_filter) & file.info(rev_filter)\$size > 0]

        if (length(fwd_ok) == 0 || length(rev_ok) == 0) {
            message("No reads passed filtering for ", sample_id, " — writing empty outputs.")
            seqtab <- matrix(integer(0), nrow=1, ncol=0, dimnames=list(sample_id, character(0)))
            read_stats <- data.frame(
                sample=sample_id, input=sum(out_filter[,1]), filtered=0L,
                denoised_fwd=0L, denoised_rev=0L, merged=0L, length_filt=0L,
                stringsAsFactors=FALSE
            )
        } else {
            # Dereplicate
            derepF <- derepFastq(fwd_ok)
            derepR <- derepFastq(rev_ok)
            if (!is.list(derepF)) derepF <- list(derepF)
            if (!is.list(derepR)) derepR <- list(derepR)

            # Sample inference
            dadaF <- dada(derepF, err = err_fwd, pool = pool, multithread = ${task.cpus})
            dadaR <- dada(derepR, err = err_rev, pool = pool, multithread = ${task.cpus})
            if (inherits(dadaF, "dada")) dadaF <- list(dadaF)
            if (inherits(dadaR, "dada")) dadaR <- list(dadaR)

            # Merge paired reads
            mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)
            if (is.data.frame(mergers)) mergers <- list(mergers)
            if (is.null(names(mergers))) names(mergers) <- sample_id

            # Make sequence table
            seqtab <- makeSequenceTable(mergers)

            # Filter by amplicon length
            seq_lengths <- nchar(getSequences(seqtab))
            seqtab <- seqtab[, seq_lengths >= min_length & seq_lengths <= max_length, drop=FALSE]

            # Read tracking stats
            read_stats <- data.frame(
                sample       = sample_id,
                input        = sum(out_filter[,1]),
                filtered     = sum(out_filter[,2]),
                denoised_fwd = sapply(dadaF, function(x) sum(x\$denoised)),
                denoised_rev = sapply(dadaR, function(x) sum(x\$denoised)),
                merged       = sapply(mergers, function(x) sum(x\$accept)),
                length_filt  = rowSums(seqtab),
                stringsAsFactors = FALSE
            )
        }

    } else {
        fwd_filter <- gsub("trimmed", "filtered", all_files)
        # See paired-end branch above for why this isn't min_length.
        min_per_read <- ifelse(trunc_len_f > 0, trunc_len_f, 20L)

        out_filter <- filterAndTrim(
            all_files, fwd_filter,
            maxEE       = max_ee_f,
            truncLen    = ifelse(trunc_len_f > 0, trunc_len_f, 0),
            truncQ      = 2,
            minLen      = min_per_read,
            rm.phix     = TRUE,
            compress    = TRUE,
            multithread = ${task.cpus}
        )

        fwd_ok <- fwd_filter[file.exists(fwd_filter) & file.info(fwd_filter)\$size > 0]

        if (length(fwd_ok) == 0) {
            message("No reads passed filtering for ", sample_id, " — writing empty outputs.")
            seqtab <- matrix(integer(0), nrow=1, ncol=0, dimnames=list(sample_id, character(0)))
            read_stats <- data.frame(
                sample=sample_id, input=sum(out_filter[,1]), filtered=0L,
                denoised=0L, length_filt=0L
            )
        } else {
            derepF <- derepFastq(fwd_ok)
            if (!is.list(derepF)) derepF <- list(derepF)

            dadaF  <- dada(derepF, err = err_fwd, pool = pool, multithread = ${task.cpus})
            if (inherits(dadaF, "dada")) dadaF <- list(dadaF)
            if (is.null(names(dadaF))) names(dadaF) <- sample_id
            seqtab <- makeSequenceTable(dadaF)

            seq_lengths <- nchar(getSequences(seqtab))
            seqtab <- seqtab[, seq_lengths >= min_length & seq_lengths <= max_length, drop=FALSE]

            read_stats <- data.frame(
                sample      = sample_id,
                input       = sum(out_filter[,1]),
                filtered    = sum(out_filter[,2]),
                denoised    = sapply(dadaF, function(x) sum(x\$denoised)),
                length_filt = rowSums(seqtab)
            )
        }
    }

    # Export
    saveRDS(seqtab, "${prefix}.asv_table.rds")

    # Write FASTA of unique ASVs
    asv_seqs   <- colnames(seqtab)
    asv_ids    <- paste0("ASV", seq_along(asv_seqs))
    fasta_out  <- paste0(">", asv_ids, "\\n", asv_seqs, collapse = "\\n")
    writeLines(fasta_out, "${prefix}.asv_seqs.fasta")

    write.table(read_stats, "${prefix}.read_stats.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)

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
    def prefix = "${meta.id}_${meta.marker}"
    """
    touch ${prefix}.asv_table.rds
    touch ${prefix}.asv_seqs.fasta
    touch ${prefix}.read_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dada2: 1.30.0
    END_VERSIONS
    """

}
