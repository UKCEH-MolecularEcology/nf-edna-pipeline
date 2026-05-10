process ECOLOGY_ALPHA {
    tag "alpha_${marker}"
    label 'process_low'

    container 'rocker/verse:4.3.3'

    publishDir "${params.outdir}/full_ecology/${marker}/01_alpha_diversity", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata
    val group_var

    output:
    tuple val(marker), path("${marker}.alpha_results/"), emit: results
    path 'versions.yml',                                  emit: versions

    script:
    def meta_arg  = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    def group_arg = group_var ? "\"${group_var}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    # ── Package loading ──────────────────────────────────────────────────────
    r_lib <- "${params.r_lib_cache}"
    dir.create(r_lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(c(r_lib, .libPaths()))
    options(mc.cores = ${task.cpus})
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


    required <- c("phyloseq","vegan","ggplot2","dplyr","tidyr","iNEXT",
                  "dunn.test","rstatix","ggpubr","cowplot","microbiome")
    invisible(lapply(required, .install_pkg))
    suppressPackageStartupMessages({
        library(phyloseq); library(vegan);  library(ggplot2)
        library(dplyr);    library(tidyr);  library(iNEXT)
        library(ggpubr);   library(cowplot); library(microbiome)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    group_var <- ${group_arg}
    out_dir   <- paste0(marker, ".alpha_results")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # ── Load data → phyloseq ─────────────────────────────────────────────
    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    tax_tab <- read.table("${taxonomy}",  sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    common  <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab <- asv_tab[common, , drop=FALSE]
    tax_tab <- tax_tab[common, , drop=FALSE]

    if (nrow(asv_tab) == 0 || ncol(asv_tab) == 0) {
        message("Empty ASV table for ", marker, " — skipping alpha diversity.")
        writeLines(c('"${task.process}":', '    skipped: empty ASV table'), "versions.yml")
        quit(status=0)
    }

    OTU <- otu_table(as.matrix(asv_tab), taxa_are_rows=TRUE)
    TAX <- tax_table(as.matrix(tax_tab))

    has_meta <- !is.null(meta_file) && file.exists(meta_file)
    if (has_meta) {
        meta    <- read.table(meta_file, sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
        # Keep only samples present in ASV table
        common_s <- intersect(rownames(meta), colnames(asv_tab))
        meta     <- meta[common_s, , drop=FALSE]
        ps       <- phyloseq(OTU, TAX, sample_data(meta))
    } else {
        ps <- phyloseq(OTU, TAX)
    }

    # Determine grouping variable
    grp <- if (!is.null(group_var) && has_meta && group_var %in% colnames(sample_data(ps))) {
        group_var
    } else if (has_meta) {
        colnames(sample_data(ps))[1]
    } else {
        NULL
    }

    # ── 1. Richness and diversity metrics ────────────────────────────────
    alpha_metrics <- c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","Fisher")
    alpha_df      <- estimate_richness(ps, measures = alpha_metrics)
    alpha_df\$sample <- rownames(alpha_df)

    # Add Pielou's evenness
    alpha_df\$Pielou_J <- alpha_df\$Shannon / log(alpha_df\$Observed)

    # Add Berger-Parker dominance
    otu_mat <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
    otu_rel  <- sweep(otu_mat, 2, colSums(otu_mat), "/")
    alpha_df\$BergerParker <- apply(otu_rel, 2, max)

    # Attach metadata
    if (has_meta) {
        meta_df <- data.frame(sample_data(ps))
        alpha_df <- merge(alpha_df, meta_df, by.x="sample", by.y="row.names", all.x=TRUE)
    }

    write.table(alpha_df, file.path(out_dir, "alpha_diversity_metrics.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # ── 2. Boxplots of alpha metrics by group ────────────────────────────
    metrics_to_plot <- c("Observed","Chao1","Shannon","InvSimpson","Pielou_J")
    plot_list <- list()

    for (m in metrics_to_plot) {
        if (!m %in% colnames(alpha_df)) next
        p <- ggplot(alpha_df, aes_string(x = if (!is.null(grp)) grp else "\"all\"",
                                          y = m,
                                          fill = if (!is.null(grp)) grp else "NULL")) +
            geom_boxplot(alpha=0.7, outlier.shape=NA) +
            geom_jitter(width=0.15, size=1.5, alpha=0.8) +
            theme_bw(base_size=12) +
            labs(title=paste(marker, "-", m), x=NULL, y=m) +
            theme(legend.position="none",
                  axis.text.x=element_text(angle=30, hjust=1))
        plot_list[[m]] <- p
    }

    combined <- plot_grid(plotlist=plot_list, ncol=3, align="hv")
    ggsave(file.path(out_dir, "alpha_diversity_boxplots.pdf"), combined, width=15, height=10)
    ggsave(file.path(out_dir, "alpha_diversity_boxplots.png"), combined, width=15, height=10, dpi=150)

    # ── 3. Statistical tests (Kruskal-Wallis + Dunn's post-hoc) ─────────
    if (!is.null(grp)) {
        stat_results <- list()
        for (m in metrics_to_plot) {
            if (!m %in% colnames(alpha_df)) next
            vals   <- alpha_df[[m]]
            groups <- alpha_df[[grp]]
            kw_res <- kruskal.test(vals ~ groups)
            stat_results[[m]] <- data.frame(
                metric    = m,
                test      = "Kruskal-Wallis",
                statistic = round(kw_res\$statistic, 4),
                df        = kw_res\$parameter,
                p_value   = round(kw_res\$p.value, 4),
                signif    = ifelse(kw_res\$p.value < 0.001, "***",
                            ifelse(kw_res\$p.value < 0.01,  "**",
                            ifelse(kw_res\$p.value < 0.05,  "*", "ns")))
            )
        }
        stat_df <- do.call(rbind, stat_results)
        write.table(stat_df, file.path(out_dir, "alpha_kruskal_wallis.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)

        # Dunn's post-hoc for significant metrics
        sig_metrics <- stat_df\$metric[stat_df\$p_value < 0.05]
        if (length(sig_metrics) > 0) {
            library(dunn.test)
            dunn_list <- list()
            for (m in sig_metrics) {
                d_res <- dunn.test(alpha_df[[m]], alpha_df[[grp]],
                                   method="bh", kw=FALSE)
                dunn_list[[m]] <- data.frame(
                    metric       = m,
                    comparison   = d_res\$comparisons,
                    Z            = round(d_res\$Z, 4),
                    p_adjusted   = round(d_res\$P.adjusted, 4)
                )
            }
            dunn_df <- do.call(rbind, dunn_list)
            write.table(dunn_df, file.path(out_dir, "alpha_dunn_posthoc.tsv"),
                        sep="\\t", quote=FALSE, row.names=FALSE)
        }
    }

    # ── 4. Rarefaction curves (iNEXT) ────────────────────────────────────
    otu_int   <- round(otu_mat)
    storage.mode(otu_int) <- "integer"
    otu_int   <- otu_int[, colSums(otu_int) > 0, drop=FALSE]

    tryCatch({
        inext_out <- iNEXT(otu_int, q=0, datatype="abundance",
                           nboot=20, conf=0.95)

        inext_df <- lapply(names(inext_out\$iNextEst\$size_based), function(sn) {
            d <- inext_out\$iNextEst\$size_based[[sn]]
            d\$sample <- sn
            d
        })
        inext_df <- do.call(rbind, inext_df)
        write.table(inext_df, file.path(out_dir, "rarefaction_curves.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)

        # Rarefaction plot
        p_rare <- ggiNEXT(inext_out, type=1) +
            theme_bw(base_size=12) +
            labs(title=paste(marker, "- Rarefaction Curves"),
                 x="Sequencing Depth", y="Species Richness") +
            theme(legend.position="right")
        ggsave(file.path(out_dir, "rarefaction_curves.pdf"), p_rare, width=10, height=7)
        ggsave(file.path(out_dir, "rarefaction_curves.png"), p_rare, width=10, height=7, dpi=150)

        # Coverage-based rarefaction
        p_cover <- ggiNEXT(inext_out, type=2) +
            theme_bw(base_size=12) +
            labs(title=paste(marker, "- Sample Coverage"),
                 x="Sequencing Depth", y="Sample Coverage") +
            theme(legend.position="right")
        ggsave(file.path(out_dir, "sample_coverage.pdf"), p_cover, width=10, height=7)
    }, error = function(e) {
        message("iNEXT failed (likely too few samples): ", conditionMessage(e))
    })

    # ── 5. Rank-abundance (Whittaker) plots ─────────────────────────────
    ra_list <- lapply(colnames(otu_rel), function(sn) {
        vals <- sort(otu_rel[, sn], decreasing=TRUE)
        vals <- vals[vals > 0]
        data.frame(sample=sn, rank=seq_along(vals), abundance=vals)
    })
    ra_df <- do.call(rbind, ra_list)
    if (has_meta) {
        meta_df  <- data.frame(sample_data(ps))
        meta_df\$sample <- rownames(meta_df)
        ra_df    <- left_join(ra_df, meta_df, by="sample")
    }
    p_ra <- ggplot(ra_df, aes(x=rank, y=log10(abundance + 1e-6),
                               group=sample,
                               color=if (!is.null(grp) && grp %in% colnames(ra_df)) ra_df[[grp]] else sample)) +
        geom_line(alpha=0.6) +
        theme_bw(base_size=12) +
        labs(title=paste(marker, "- Rank-Abundance (Whittaker)"),
             x="Rank", y="log10(Relative Abundance)", color=grp)
    ggsave(file.path(out_dir, "rank_abundance.pdf"), p_ra, width=10, height=6)
    ggsave(file.path(out_dir, "rank_abundance.png"), p_ra, width=10, height=6, dpi=150)

    # ── 6. Occupancy-abundance plot ──────────────────────────────────────
    occupancy  <- rowMeans(otu_mat > 0)
    mean_abund <- rowMeans(otu_rel)
    oa_df <- data.frame(
        asv_id    = rownames(otu_mat),
        occupancy = occupancy,
        mean_rel_abund = mean_abund
    )
    if ("Genus" %in% colnames(tax_tab)) {
        oa_df\$genus <- tax_tab[oa_df\$asv_id, "Genus"]
    }
    write.table(oa_df, file.path(out_dir, "occupancy_abundance.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    p_oa <- ggplot(oa_df, aes(x=occupancy, y=log10(mean_rel_abund + 1e-6))) +
        geom_point(alpha=0.4, size=1.5) +
        geom_smooth(method="loess", se=TRUE, color="steelblue") +
        theme_bw(base_size=12) +
        labs(title=paste(marker, "- Occupancy-Abundance"),
             x="Occupancy (fraction of samples)", y="log10(Mean Rel. Abundance)")
    ggsave(file.path(out_dir, "occupancy_abundance.pdf"), p_oa, width=8, height=6)

    message("Alpha diversity analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    phyloseq: ', packageVersion('phyloseq')),
        paste0('    iNEXT: ',    packageVersion('iNEXT')),
        paste0('    vegan: ',    packageVersion('vegan')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
