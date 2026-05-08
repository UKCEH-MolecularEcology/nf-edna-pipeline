process MERGE_ASV_TABLES {
    tag "merge_${marker}"
    label 'process_medium'

    container 'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'

    publishDir "${params.outdir}/asv_tables/${marker}", mode: 'copy'

    input:
    tuple val(marker), val(metas), path(rds_files)

    output:
    tuple val(marker), path('*.merged_asv_table.tsv'),  emit: merged_table
    tuple val(marker), path('*.merged_asv_table.rds'),  emit: merged_rds
    tuple val(marker), path('*.read_tracking.tsv'),     emit: read_tracking
    path 'versions.yml',                                 emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)

    marker  <- "${marker}"
    rds_files <- list.files(".", pattern = "\\.asv_table\\.rds\$", full.names = TRUE)

    # Merge sequence tables across samples
    seqtabs <- lapply(rds_files, readRDS)
    merged  <- mergeSequenceTables(seqtabs)

    # Remove bimeras (chimeras) from the merged table
    # (De novo chimera check on merged table as final step)
    merged_nochim <- removeBimeraDenovo(
        merged,
        method      = "${params.dada2_chimera}",
        multithread = ${task.cpus},
        verbose     = TRUE
    )

    # Read tracking: how many reads survived each step
    getN <- function(x) sum(getUniques(x))

    read_track <- data.frame(
        sample     = rownames(merged_nochim),
        merged_in  = rowSums(merged),
        nonchimera = rowSums(merged_nochim)
    )
    write.table(read_track, paste0(marker, ".read_tracking.tsv"),
                sep = "\\t", quote = FALSE, row.names = FALSE)

    # Export as TSV (samples x ASVs)
    asv_tab <- as.data.frame(t(merged_nochim))
    asv_tab <- cbind(asv_id = rownames(asv_tab), asv_tab)
    write.table(asv_tab, paste0(marker, ".merged_asv_table.tsv"),
                sep = "\\t", quote = FALSE, row.names = FALSE)

    saveRDS(merged_nochim, paste0(marker, ".merged_asv_table.rds"))

    writeLines(
        c(
            paste0('"${task.process}":'),
            paste0('    dada2: ', packageVersion('dada2')),
            paste0('    R: ', R.version\$major, '.', R.version\$minor)
        ),
        "versions.yml"
    )
    """
}
