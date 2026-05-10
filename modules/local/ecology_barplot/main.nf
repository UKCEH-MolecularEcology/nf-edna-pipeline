process ECOLOGY_BARPLOT {
    tag "barplot_${marker}"
    label 'process_low'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/ecology:r4.3.3'

    publishDir "${params.outdir}/ecology/${marker}/composition", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata

    output:
    tuple val(marker), path("${marker}.composition_results/"), emit: results
    path 'versions.yml',                                emit: versions

    script:
    def meta_arg = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    # ── Package loading ──────────────────────────────────────────────────────
    r_lib <- "${params.r_lib_cache}"
    dir.create(r_lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(c(r_lib, .libPaths()))
    options(mc.cores = ${task.cpus})
    if (length(readLines("${asv_table}")) <= 1L) {
        dir.create("${marker}.composition_results", showWarnings=FALSE, recursive=TRUE)
        writeLines("skipped: empty ASV table", "${marker}.composition_results/skipped.txt")
        writeLines(c('"${task.process}":', '    skipped: empty ASV table'), "versions.yml")
        quit(status=0)
    }
    .install_pkg <- function(pkg) {
        if (requireNamespace(pkg, quietly=TRUE)) return(invisible(NULL))
        install.packages(pkg, quiet=TRUE)
        if (!requireNamespace(pkg, quietly=TRUE)) {
            Sys.sleep(runif(1, 5, 15))
            install.packages(pkg, quiet=TRUE, INSTALL_opts="--no-lock")
        }
    }

    r_ver <- numeric_version(paste(R.version\$major, R.version\$minor, sep="."))
    bioc_ver <- if (r_ver >= "4.5") "3.22" else if (r_ver >= "4.4") "3.20" else if (r_ver >= "4.3") "3.18" else "3.16"
    options(repos = c(
        BioCsoft = paste0("https://bioconductor.org/packages/", bioc_ver, "/bioc"),
        BioCann  = paste0("https://bioconductor.org/packages/", bioc_ver, "/data/annotation"),
        CRAN     = "https://cloud.r-project.org"
    ))


    required <- c("phyloseq", "ggplot2", "dplyr", "tidyr")
    invisible(lapply(required, .install_pkg))
    suppressPackageStartupMessages({
        library(phyloseq)
        library(ggplot2)
        library(dplyr)
        library(tidyr)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    out_dir   <- paste0(marker, ".composition_results")
    dir.create(out_dir, showWarnings = FALSE)

    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE,
                          row.names=1, check.names=FALSE)
    asv_tab[is.na(asv_tab)] <- 0
    tax_tab <- read.table("${taxonomy}", sep="\\t", header=TRUE,
                          row.names=1, check.names=FALSE)

    common_asvs <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab     <- asv_tab[common_asvs, , drop=FALSE]
    tax_tab     <- tax_tab[common_asvs, , drop=FALSE]

    if (nrow(asv_tab) == 0 || ncol(asv_tab) == 0) {
        message("Empty ASV table for ", marker, " — skipping composition analysis.")
        writeLines(c('"${task.process}":', '    skipped: empty ASV table'), "versions.yml")
        writeLines("skipped: empty ASV table", file.path(out_dir, "skipped.txt"))
        quit(status=0)
    }

    OTU <- otu_table(as.matrix(asv_tab), taxa_are_rows = TRUE)
    TAX <- tax_table(as.matrix(tax_tab))

    if (!is.null(meta_file) && file.exists(meta_file)) {
        meta  <- read.table(meta_file, sep="\\t", header=TRUE,
                            row.names=1, check.names=FALSE)
        ps    <- phyloseq(OTU, TAX, sample_data(meta))
    } else {
        ps    <- phyloseq(OTU, TAX)
    }

    ps <- prune_samples(sample_sums(ps) > 0, ps)
    ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    if (nsamples(ps) == 0 || ntaxa(ps) == 0) {
        writeLines("skipped: no data after filtering", file.path(out_dir, "skipped.txt"))
        writeLines(c(paste0('"${task.process}":'), '    skipped: no data after filtering'), "versions.yml")
        quit(status=0)
    }

    # Transform to relative abundance
    ps_rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)

    # Determine rank levels based on what's available
    available_ranks <- c("Phylum","Class","Order","Family","Genus")
    tax_ranks <- intersect(available_ranks, colnames(tax_table(ps)))

    for (rank in tax_ranks) {
        tryCatch({
        ps_rank <- tax_glom(ps_rel, taxrank = rank, NArm = FALSE)
        if (ntaxa(ps_rank) == 0) stop("no taxa")

        # Keep top N taxa, lump rest as "Other"
        top_n   <- 20
        top_ids <- names(sort(taxa_sums(ps_rank), decreasing=TRUE))[seq_len(min(top_n, ntaxa(ps_rank)))]
        ps_top  <- prune_taxa(top_ids, ps_rank)

        # Melt to long format
        df <- psmelt(ps_top)
        df[[rank]][is.na(df[[rank]])] <- "Unclassified"

        # Stacked barplot
        p <- ggplot(df, aes(x=Sample, y=Abundance, fill=.data[[rank]])) +
            geom_bar(stat="identity") +
            scale_y_continuous(labels = function(x) paste0(x, "%")) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle=45, hjust=1, size=8),
                legend.text = element_text(size=7),
                legend.key.size = unit(0.4, "cm")
            ) +
            labs(title = paste(marker, "-", rank, "composition"),
                 y = "Relative Abundance (%)",
                 x = "Sample")

        ggsave(file.path(out_dir, paste0(rank, "_barplot.pdf")), p,
               width=max(8, nsamples(ps)*0.4 + 3), height=7)
        ggsave(file.path(out_dir, paste0(rank, "_barplot.png")), p,
               width=max(8, nsamples(ps)*0.4 + 3), height=7, dpi=150)

        # Export composition table
        comp_wide <- df %>%
            select(Sample, !!sym(rank), Abundance) %>%
            pivot_wider(names_from=Sample, values_from=Abundance, values_fill=0)
        write.table(comp_wide,
                    file.path(out_dir, paste0(rank, "_relative_abundance.tsv")),
                    sep="\\t", quote=FALSE, row.names=FALSE)
        }, error=function(e) message("Rank ", rank, " skipped: ", conditionMessage(e)))
    }

    # Heatmap of top ASVs
    if (ntaxa(ps_rel) > 0) {
    top_asvs <- names(sort(taxa_sums(ps_rel), decreasing=TRUE))[seq_len(min(50, ntaxa(ps_rel)))]
    ps_heat  <- prune_taxa(top_asvs, ps_rel)

    heat_mat <- as.matrix(t(otu_table(ps_heat)))
    tax_labels <- if ("Genus" %in% colnames(tax_table(ps_heat))) {
        paste0(tax_table(ps_heat)[,"Genus"], "_", rownames(tax_table(ps_heat)))
    } else {
        rownames(tax_table(ps_heat))
    }
    colnames(heat_mat) <- tax_labels

    if (nrow(heat_mat) >= 2 && ncol(heat_mat) >= 2) {
    pdf(file.path(out_dir, "top_asvs_heatmap.pdf"), width=14, height=8)
    heatmap(heat_mat, Rowv=NA, Colv=NA, col=heat.colors(256),
            margins=c(12,6), main=paste(marker, "- Top ASVs (Rel. Abund. %)"),
            cexCol=0.6, cexRow=0.8)
    dev.off()
    }

    } # end ntaxa > 0
    message("Composition analysis complete.")

    writeLines(
        c(
            paste0('"${task.process}":'),
            paste0('    phyloseq: ', packageVersion('phyloseq')),
            paste0('    ggplot2: ',  packageVersion('ggplot2')),
            paste0('    R: ', R.version\$major, '.', R.version\$minor)
        ),
        "versions.yml"
    )
    """
}
