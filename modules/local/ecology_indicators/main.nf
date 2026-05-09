process ECOLOGY_INDICATORS {
    tag "indicators_${marker}"
    label 'process_medium'

    container 'ghcr.io/rocker-project/verse:4.3.3'

    publishDir "${params.outdir}/full_ecology/${marker}/06_indicator_species", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata
    val group_var

    output:
    tuple val(marker), path("${marker}.indicator_results/"), emit: results
    path 'versions.yml',                                      emit: versions

    script:
    def meta_arg  = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    def group_arg = group_var ? "\"${group_var}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    pkgs <- c("vegan","indicspecies","ggplot2","dplyr","tidyr","microbiome","phyloseq")
    for (pkg in pkgs) {
        if (!requireNamespace(pkg, quietly=TRUE)) {
            if (pkg %in% c("microbiome","phyloseq")) {
                if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
                BiocManager::install(pkg, update=FALSE, ask=FALSE)
            } else {
                install.packages(pkg, repos="https://cloud.r-project.org")
            }
        }
    }
    suppressPackageStartupMessages({
        library(vegan); library(indicspecies); library(ggplot2)
        library(dplyr); library(tidyr);        library(microbiome)
        library(phyloseq)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    group_var <- ${group_arg}
    out_dir   <- paste0(marker, ".indicator_results")
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    set.seed(42)

    # ── Load data ────────────────────────────────────────────────────────
    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    tax_tab <- read.table("${taxonomy}",  sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    common  <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab <- asv_tab[common, , drop=FALSE]
    tax_tab <- tax_tab[common, , drop=FALSE]

    has_meta <- !is.null(meta_file) && file.exists(meta_file)
    meta     <- if (has_meta) {
        read.table(meta_file, sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    } else {
        NULL
    }

    grp <- if (!is.null(meta) && !is.null(group_var) && group_var %in% colnames(meta)) {
        group_var
    } else if (!is.null(meta)) {
        cat_cols <- names(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
        if (length(cat_cols) > 0) cat_cols[1] else NULL
    } else {
        NULL
    }

    otu_t <- t(asv_tab)

    # Tax label helper
    get_label <- function(asv_ids) {
        ranks <- c("Genus","Family","Order","Class","Phylum")
        avail <- intersect(ranks, colnames(tax_tab))
        if (length(avail) == 0) return(asv_ids)
        tax_tab[asv_ids, avail[1]]
    }

    # ── 1. IndVal — Indicator Species Analysis ────────────────────────────
    if (!is.null(grp)) {
        common_s <- intersect(rownames(meta), rownames(otu_t))
        groups   <- as.character(meta[common_s, grp])
        otu_sub  <- otu_t[common_s, , drop=FALSE]

        tryCatch({
            indval_res <- multipatt(otu_sub, groups,
                                    control=how(nperm=999),
                                    func="IndVal.g", duleg=FALSE)
            sink(file.path(out_dir, "indval_summary.txt"))
            summary(indval_res, alpha=0.05, indvalcomp=TRUE)
            sink()

            # Format as data frame
            iv_df <- data.frame(
                asv_id   = rownames(indval_res\$sign),
                indval_res\$sign,
                stringsAsFactors = FALSE
            )
            iv_df\$tax_label <- get_label(iv_df\$asv_id)
            iv_df <- iv_df[!is.na(iv_df\$p.value) & iv_df\$p.value < 0.05, ]
            iv_df <- iv_df[order(iv_df\$p.value), ]
            write.table(iv_df, file.path(out_dir, "indval_significant.tsv"),
                        sep="\\t", quote=FALSE, row.names=FALSE)

            # Plot top indicator species
            if (nrow(iv_df) > 0) {
                top_iv <- head(iv_df, 30)
                p_iv <- ggplot(top_iv, aes(x=reorder(tax_label, stat),
                                            y=stat, fill=p.value < 0.05)) +
                    geom_col() +
                    coord_flip() +
                    scale_fill_manual(values=c("FALSE"="grey70","TRUE"="#27AE60")) +
                    theme_bw(base_size=11) +
                    labs(title=paste(marker, "- IndVal Indicator Species (p < 0.05)"),
                         x=NULL, y="IndVal statistic", fill="Significant")
                ggsave(file.path(out_dir, "indval_barplot.pdf"), p_iv,
                       width=10, height=max(5, nrow(top_iv)*0.3))
                ggsave(file.path(out_dir, "indval_barplot.png"), p_iv,
                       width=10, height=max(5, nrow(top_iv)*0.3), dpi=150)
            }
        }, error=function(e) message("IndVal failed: ", conditionMessage(e)))

        # ── 2. SIMPER ──────────────────────────────────────────────────────
        tryCatch({
            simp_res <- simper(otu_sub, groups, permutations=999)
            sink(file.path(out_dir, "simper_summary.txt"))
            print(summary(simp_res, ordered=TRUE))
            sink()

            # Export per-comparison SIMPER tables
            for (comp_name in names(simp_res)) {
                s <- summary(simp_res)[[comp_name]]
                if (!is.null(s) && nrow(s) > 0) {
                    s\$asv_id    <- rownames(s)
                    s\$tax_label <- get_label(s\$asv_id)
                    s\$comparison <- comp_name
                    write.table(s, file.path(out_dir,
                        paste0("simper_", gsub(" ", "_", comp_name), ".tsv")),
                        sep="\\t", quote=FALSE, row.names=FALSE, na="")
                }
            }

            # Combined plot: top contributors for first comparison
            first_comp <- names(simp_res)[1]
            s_first    <- summary(simp_res)[[first_comp]]
            if (!is.null(s_first) && nrow(s_first) > 0) {
                s_first\$asv_id    <- rownames(s_first)
                s_first\$tax_label <- get_label(s_first\$asv_id)
                top_simp <- head(s_first[order(-s_first\$average), ], 20)
                p_simp <- ggplot(top_simp, aes(x=reorder(tax_label, average), y=average*100)) +
                    geom_col(fill="#2980B9", alpha=0.85) +
                    geom_col(aes(y=cumsum(average*100)), fill=NA, color="black", linetype=2, linewidth=0.3) +
                    coord_flip() +
                    theme_bw(base_size=11) +
                    labs(title=paste(marker, "- SIMPER:", first_comp),
                         x=NULL, y="Average contribution (%)")
                ggsave(file.path(out_dir, paste0("simper_", gsub(" ","_",first_comp), "_plot.pdf")),
                       p_simp, width=10, height=max(5, nrow(top_simp)*0.35))
            }
        }, error=function(e) message("SIMPER failed: ", conditionMessage(e)))
    }

    # ── 3. Core microbiome ────────────────────────────────────────────────
    # Core = ASVs present in >= 50% of samples at >= 0.1% relative abundance
    otu_rel <- sweep(asv_tab, 2, colSums(asv_tab), "/") * 100
    core_matrix <- data.frame(
        asv_id        = rownames(otu_rel),
        mean_rel_abund = rowMeans(otu_rel),
        prevalence     = rowMeans(asv_tab > 0) * 100
    )
    core_matrix\$tax_label <- get_label(core_matrix\$asv_id)

    for (det_thresh in c(0.01, 0.1, 1)) {
        for (prev_thresh in c(50, 75, 90)) {
            is_core <- core_matrix\$mean_rel_abund >= det_thresh &
                       core_matrix\$prevalence      >= prev_thresh
            core_matrix[[paste0("core_", det_thresh, "_pct_det_", prev_thresh, "_pct_prev")]] <- is_core
        }
    }
    write.table(core_matrix[order(-core_matrix\$prevalence), ],
                file.path(out_dir, "core_microbiome_matrix.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # Core occupancy-abundance plot
    p_core <- ggplot(core_matrix, aes(x=prevalence, y=log10(mean_rel_abund + 0.001))) +
        geom_point(aes(color=core_matrix[[paste0("core_0.1_pct_det_50_pct_prev")]]),
                   alpha=0.5, size=1.5) +
        scale_color_manual(values=c("FALSE"="grey70","TRUE"="#E74C3C")) +
        geom_vline(xintercept=50, linetype="dashed", color="grey40") +
        geom_hline(yintercept=log10(0.1), linetype="dashed", color="grey40") +
        theme_bw(base_size=12) +
        labs(title=paste(marker, "- Core Microbiome"),
             subtitle="Red = core (prevalence ≥ 50%, det ≥ 0.1%)",
             x="Prevalence (%)", y="log10(Mean Relative Abundance %)",
             color="Core")
    ggsave(file.path(out_dir, "core_microbiome_plot.pdf"), p_core, width=8, height=6)
    ggsave(file.path(out_dir, "core_microbiome_plot.png"), p_core, width=8, height=6, dpi=150)

    # ── 4. Taxon prevalence × mean abundance heatmap ──────────────────────
    if (!is.null(grp) && !is.null(meta)) {
        common_s <- intersect(rownames(meta), colnames(asv_tab))
        groups   <- as.character(meta[common_s, grp])
        otu_sub  <- asv_tab[, common_s, drop=FALSE]
        otu_rel_sub <- sweep(otu_sub, 2, colSums(otu_sub), "/") * 100

        group_stats <- lapply(unique(groups), function(g) {
            idx   <- which(groups == g)
            m_rel <- rowMeans(otu_rel_sub[, idx, drop=FALSE])
            prev  <- rowMeans(otu_sub[, idx, drop=FALSE] > 0) * 100
            data.frame(asv_id=rownames(otu_rel_sub), group=g,
                       mean_rel=m_rel, prevalence=prev)
        })
        group_df <- do.call(rbind, group_stats)
        write.table(group_df, file.path(out_dir, "per_group_taxon_stats.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
    }

    message("Indicator species analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    indicspecies: ', packageVersion('indicspecies')),
        paste0('    vegan: ',        packageVersion('vegan')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
