process DADA2_DENOISE {
    tag "${meta.id}_${marker}"
    label 'process_high'

    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/dada2/${marker}/asv_tables", mode: 'copy'

    input:
    tuple val(run_id), val(meta), path(reads), path(err_fwd), path(err_rev), path(all_filtered), path(filter_stats)
    val marker
    val pool
    val min_length
    val max_length

    output:
    tuple val(meta), path('*.asv_table.rds'),   emit: asv_table
    tuple val(meta), path('*.asv_seqs.fasta'),  emit: asv_seqs
    tuple val(meta), path('*.read_stats.tsv'),  emit: read_stats
    path 'versions.yml',                         emit: versions

    script:
    def prefix   = "${meta.id}_${marker}"
    def pe       = meta.single_end ? 'FALSE' : 'TRUE'
    def pool_arg = pool == true ? 'TRUE' : pool == 'pseudo' ? '"pseudo"' : 'FALSE'
    """
    #!/usr/bin/env Rscript
    library(dada2)

    sample_id   <- "${meta.id}"
    marker      <- "${marker}"
    paired_end  <- as.logical("${pe}")
    pool        <- ${pool_arg}
    min_length  <- as.integer("${min_length}")
    max_length  <- as.integer("${max_length}")

    err_fwd <- readRDS("${err_fwd}")
    err_rev <- readRDS("${err_rev}")

    # Reuse the files DADA2_LEARN_ERRORS already filtered for the whole run
    # (filtering every sample twice -- once collectively there, once again
    # per-sample here -- doubled the total filtering work for no benefit).
    filter_stats <- read.table("${filter_stats}", sep = "\\t", header = TRUE, check.names = FALSE)
    stats_row <- function(fname) {
        r <- filter_stats[filter_stats\$file == fname, , drop = FALSE]
        if (nrow(r) == 1) r else NULL
    }

    if (paired_end) {
        fwd_ok <- paste0(sample_id, "_", marker, "_R1.filtered.fastq.gz")
        rev_ok <- paste0(sample_id, "_", marker, "_R2.filtered.fastq.gz")
        fwd_row <- stats_row(sub("filtered", "trimmed", fwd_ok))
        n_input <- if (!is.null(fwd_row)) fwd_row\$reads.in else NA_integer_

        if (!file.exists(fwd_ok) || !file.exists(rev_ok) ||
            file.info(fwd_ok)\$size == 0 || file.info(rev_ok)\$size == 0) {
            message("No reads passed filtering for ", sample_id, " — writing empty outputs.")
            seqtab <- matrix(integer(0), nrow=1, ncol=0, dimnames=list(sample_id, character(0)))
            read_stats <- data.frame(
                sample=sample_id, input=n_input, filtered=0L,
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
                input        = n_input,
                filtered     = if (!is.null(fwd_row)) fwd_row\$reads.out else sum(sapply(derepF, function(x) sum(x\$uniques))),
                denoised_fwd = sapply(dadaF, function(x) sum(x\$denoised)),
                denoised_rev = sapply(dadaR, function(x) sum(x\$denoised)),
                merged       = sapply(mergers, function(x) sum(x\$accept)),
                length_filt  = rowSums(seqtab),
                stringsAsFactors = FALSE
            )
        }

    } else {
        fwd_ok  <- paste0(sample_id, "_", marker, "_R1.filtered.fastq.gz")
        fwd_row <- stats_row(sub("filtered", "trimmed", fwd_ok))
        n_input <- if (!is.null(fwd_row)) fwd_row\$reads.in else NA_integer_

        if (!file.exists(fwd_ok) || file.info(fwd_ok)\$size == 0) {
            message("No reads passed filtering for ", sample_id, " — writing empty outputs.")
            seqtab <- matrix(integer(0), nrow=1, ncol=0, dimnames=list(sample_id, character(0)))
            read_stats <- data.frame(
                sample=sample_id, input=n_input, filtered=0L,
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
                input       = n_input,
                filtered    = if (!is.null(fwd_row)) fwd_row\$reads.out else sum(sapply(derepF, function(x) sum(x\$uniques))),
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
