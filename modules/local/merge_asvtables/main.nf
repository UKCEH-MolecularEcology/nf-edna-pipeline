process MERGE_ASV_TABLES {
    tag "merge_${marker}"
    label 'process_medium'

    container 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/asv_tables/${marker}", mode: 'copy'

    input:
    tuple val(marker), val(metas), path(rds_files)

    output:
    tuple val(marker), path('*.merged_asv_table.tsv'),  emit: merged_table
    tuple val(marker), path('*.merged_asv_table.rds'),  emit: merged_rds
    tuple val(marker), path('*.read_tracking.tsv'),     emit: read_tracking
    tuple val(marker), path('*.merged_asv.fasta'),      emit: merged_fasta
    tuple val(marker), path('*.asv_lookup.tsv'),        emit: asv_lookup
    path 'versions.yml',                                 emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)

    marker  <- "${marker}"
    rds_files <- list.files(".", pattern = "\\\\.asv_table\\\\.rds\$", full.names = TRUE)

    # Merge sequence tables across samples (skip empty tables from zero-read samples)
    seqtabs <- lapply(rds_files, readRDS)
    seqtabs <- Filter(function(x) is.matrix(x) && ncol(x) > 0 && nrow(x) > 0, seqtabs)

    if (length(seqtabs) == 0) {
        message("No non-empty sequence tables found for ", marker, " — writing empty outputs.")
        empty_mat <- matrix(integer(0), nrow=0, ncol=0)
        saveRDS(empty_mat, paste0(marker, ".merged_asv_table.rds"))
        write.table(data.frame(asv_id=character(0)), paste0(marker, ".merged_asv_table.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
        write.table(data.frame(sample=character(0), merged_in=integer(0), nonchimera=integer(0)),
                    paste0(marker, ".read_tracking.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
        file.create(paste0(marker, ".merged_asv.fasta"))
        write.table(data.frame(ASV=character(0), sequence=character(0), length=integer(0)),
                    paste0(marker, ".asv_lookup.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
    } else {
        merged <- if (length(seqtabs) == 1) seqtabs[[1]] else mergeSequenceTables(tables = seqtabs)

        merged_nochim <- removeBimeraDenovo(
            merged,
            method      = "${params.dada2_chimera}",
            multithread = ${task.cpus},
            verbose     = TRUE
        )

        read_track <- data.frame(
            sample     = rownames(merged_nochim),
            merged_in  = rowSums(merged),
            nonchimera = rowSums(merged_nochim)
        )
        write.table(read_track, paste0(marker, ".read_tracking.tsv"),
                    sep = "\\t", quote = FALSE, row.names = FALSE)

        asv_tab <- as.data.frame(t(merged_nochim))
        asv_tab <- cbind(asv_id = rownames(asv_tab), asv_tab)
        write.table(asv_tab, paste0(marker, ".merged_asv_table.tsv"),
                    sep = "\\t", quote = FALSE, row.names = FALSE)

        saveRDS(merged_nochim, paste0(marker, ".merged_asv_table.rds"))

        # ASV FASTA + ID<->sequence lookup, keyed by sequential ASV IDs in
        # column order — needed by marker-specific downstream taxonomy
        # branches (e.g. 12S SINTAX) that require a query FASTA rather than
        # a sequence-keyed table.
        asv_seqs <- colnames(merged_nochim)
        asv_ids  <- paste0("ASV", seq_along(asv_seqs))

        fasta_lines <- character(2 * length(asv_seqs))
        fasta_lines[c(TRUE, FALSE)] <- paste0(">", asv_ids)
        fasta_lines[c(FALSE, TRUE)] <- asv_seqs
        writeLines(fasta_lines, paste0(marker, ".merged_asv.fasta"))

        write.table(
            data.frame(ASV = asv_ids, sequence = asv_seqs, length = nchar(asv_seqs)),
            paste0(marker, ".asv_lookup.tsv"),
            sep = "\\t", quote = FALSE, row.names = FALSE
        )
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
    touch ${marker}.merged_asv_table.tsv
    touch ${marker}.merged_asv_table.rds
    touch ${marker}.read_tracking.tsv
    touch ${marker}.merged_asv.fasta
    touch ${marker}.asv_lookup.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.3
    END_VERSIONS
    """
}
