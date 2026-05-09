process ECOLOGY_DIFFERENTIAL {
    tag "diffabund_${marker}"
    label 'process_high'

    container 'bioconductor/bioconductor_docker:RELEASE_3_18'

    publishDir "${params.outdir}/full_ecology/${marker}/04_differential_abundance", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata
    val group_var

    output:
    tuple val(marker), path("${marker}.diffabund_results/"), emit: results
    path 'versions.yml',                                      emit: versions

    when:
    metadata.name != 'NO_FILE'

    script:
    def meta_arg  = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    def group_arg = group_var ? "\"${group_var}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    r_lib <- file.path(getwd(), ".r_libs")
    dir.create(r_lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(c(r_lib, .libPaths()))
    install.packages("BiocManager", repos="https://cloud.r-project.org",
                     quiet=TRUE, lib=r_lib)

    pkgs <- c("DESeq2","ALDEx2","ggplot2","dplyr","tidyr","phyloseq",
              "ggrepel","ComplexHeatmap")
    for (pkg in pkgs) {
        if (!requireNamespace(pkg, quietly=TRUE)) {
            BiocManager::install(pkg, update=FALSE, ask=FALSE)
        }
    }
    suppressPackageStartupMessages({
        library(DESeq2); library(ALDEx2); library(ggplot2)
        library(dplyr);  library(phyloseq); library(ggrepel)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    group_var <- ${group_arg}
    out_dir   <- paste0(marker, ".diffabund_results")
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    set.seed(42)

    if (is.null(meta_file) || !file.exists(meta_file)) {
        message("No metadata provided — skipping differential abundance analysis.")
        writeLines(c('"${task.process}":', '    skipped: no metadata'), "versions.yml")
        quit(status=0)
    }

    # ── Load data ────────────────────────────────────────────────────────
    asv_tab  <- read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    tax_tab  <- read.table("${taxonomy}",  sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    meta     <- read.table(meta_file, sep="\\t", header=TRUE, row.names=1, check.names=FALSE)

    common_asvs <- intersect(rownames(asv_tab), rownames(tax_tab))
    common_samp <- intersect(colnames(asv_tab), rownames(meta))
    asv_tab  <- asv_tab[common_asvs, common_samp, drop=FALSE]
    tax_tab  <- tax_tab[common_asvs, , drop=FALSE]
    meta     <- meta[common_samp, , drop=FALSE]

    # Determine grouping variable
    grp <- if (!is.null(group_var) && group_var %in% colnames(meta)) {
        group_var
    } else {
        # Use first categorical column
        cat_cols <- names(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
        if (length(cat_cols) > 0) cat_cols[1] else NULL
    }

    if (is.null(grp)) {
        message("No categorical grouping variable found. Skipping differential abundance.")
        writeLines(c('"${task.process}":', '    skipped: no group variable'), "versions.yml")
        quit(status=0)
    }

    meta[[grp]] <- factor(meta[[grp]])
    grp_levels  <- levels(meta[[grp]])

    if (length(grp_levels) < 2) {
        message("Need at least 2 groups. Found: ", length(grp_levels))
        writeLines(c('"${task.process}":', '    skipped: <2 groups'), "versions.yml")
        quit(status=0)
    }

    # Remove ASVs present in <10% of samples (reduce sparsity)
    min_prev <- ceiling(0.1 * ncol(asv_tab))
    asv_tab  <- asv_tab[rowSums(asv_tab > 0) >= min_prev, , drop=FALSE]

    # ── Helper: add taxonomy labels ─────────────────────────────────────
    add_tax_label <- function(df) {
        tax_ranks <- c("Genus","Family","Order","Class","Phylum")
        avail     <- intersect(tax_ranks, colnames(tax_tab))
        if (length(avail) == 0) {
            df\$tax_label <- df\$asv_id
            return(df)
        }
        top_rank <- avail[1]
        labels   <- tax_tab[df\$asv_id, top_rank]
        labels[is.na(labels) | labels == ""] <- df\$asv_id[is.na(labels) | labels == ""]
        df\$tax_label <- labels
        df
    }

    # ── Helper: volcano plot ─────────────────────────────────────────────
    make_volcano <- function(df, title, fc_col, pval_col) {
        df\$significant <- df[[pval_col]] < 0.05 & abs(df[[fc_col]]) > 1
        df\$direction   <- ifelse(df[[fc_col]] > 1, "Up", ifelse(df[[fc_col]] < -1, "Down", "NS"))
        df\$direction[!df\$significant] <- "NS"

        ggplot(df, aes_string(x=fc_col, y=paste0("-log10(", pval_col, " + 1e-300)"))) +
            geom_point(aes(color=direction), alpha=0.5, size=1.5) +
            scale_color_manual(values=c("Up"="#C0392B","Down"="#2980B9","NS"="grey60")) +
            geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey40") +
            geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
            geom_text_repel(
                data=subset(df, significant),
                aes_string(label="tax_label"), size=2.5, max.overlaps=15
            ) +
            theme_bw(base_size=12) +
            labs(title=title, x="log2 Fold Change",
                 y="-log10(adj. p-value)", color="Direction")
    }

    all_pairwise <- combn(grp_levels, 2, simplify=FALSE)

    # ── DESeq2 ──────────────────────────────────────────────────────────
    message("Running DESeq2...")
    asv_int <- round(asv_tab)
    storage.mode(asv_int) <- "integer"

    dds <- DESeqDataSetFromMatrix(
        countData = asv_int,
        colData   = meta[colnames(asv_int), , drop=FALSE],
        design    = as.formula(paste0("~ ", grp))
    )
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds, quiet=TRUE, fitType="local")

    deseq_all <- list()
    for (pair in all_pairwise) {
        contrast <- c(grp, pair[2], pair[1])
        tryCatch({
            res <- results(dds, contrast=contrast, alpha=0.05)
            res_df <- data.frame(
                asv_id  = rownames(res),
                log2FC  = round(res\$log2FoldChange, 4),
                baseMean = round(res\$baseMean, 2),
                lfcSE   = round(res\$lfcSE, 4),
                stat    = round(res\$stat, 4),
                pvalue  = res\$pvalue,
                padj    = res\$padj
            )
            res_df <- add_tax_label(res_df)
            res_df\$comparison <- paste0(pair[2], "_vs_", pair[1])

            deseq_all[[paste0(pair[2], "_vs_", pair[1])]] <- res_df

            out_file <- file.path(out_dir, paste0("deseq2_", pair[2], "_vs_", pair[1], ".tsv"))
            write.table(res_df[order(res_df\$padj), ],
                        out_file, sep="\\t", quote=FALSE, row.names=FALSE, na="")

            # Volcano plot
            res_df_nona <- res_df[!is.na(res_df\$padj), ]
            if (nrow(res_df_nona) > 0) {
                p_vol <- make_volcano(res_df_nona,
                    paste("DESeq2:", pair[2], "vs", pair[1], "(", marker, ")"),
                    "log2FC", "padj")
                ggsave(file.path(out_dir, paste0("volcano_deseq2_", pair[2], "_vs_", pair[1], ".pdf")),
                       p_vol, width=9, height=7)
                ggsave(file.path(out_dir, paste0("volcano_deseq2_", pair[2], "_vs_", pair[1], ".png")),
                       p_vol, width=9, height=7, dpi=150)
            }
        }, error=function(e) message("DESeq2 failed for ", pair[2], " vs ", pair[1], ": ", conditionMessage(e)))
    }

    # MA plot for first comparison
    if (length(deseq_all) > 0) {
        first_res <- deseq_all[[1]]
        first_res\$sig <- !is.na(first_res\$padj) & first_res\$padj < 0.05
        p_ma <- ggplot(first_res, aes(x=log2(baseMean+1), y=log2FC, color=sig)) +
            geom_point(alpha=0.4, size=1) +
            scale_color_manual(values=c("FALSE"="grey60","TRUE"="#C0392B")) +
            geom_hline(yintercept=0, linetype="dashed") +
            theme_bw(base_size=12) +
            labs(title=paste("DESeq2 MA Plot:", names(deseq_all)[1], "-", marker),
                 x="log2(Mean Counts + 1)", y="log2 Fold Change", color="Significant")
        ggsave(file.path(out_dir, "ma_plot_deseq2.pdf"), p_ma, width=8, height=6)
    }

    # ── ALDEx2 (compositional approach) ─────────────────────────────────
    message("Running ALDEx2...")
    if (length(grp_levels) == 2) {
        tryCatch({
            conds <- as.character(meta[colnames(asv_int), grp])
            aldex_out <- aldex(asv_int, conds, mc.samples=128,
                               test="t", effect=TRUE, denom="all")
            aldex_df <- data.frame(
                asv_id    = rownames(aldex_out),
                effect    = round(aldex_out\$effect, 4),
                overlap   = round(aldex_out\$overlap, 4),
                we.ep     = aldex_out\$we.ep,
                we.eBH    = aldex_out\$we.eBH,
                wi.ep     = aldex_out\$wi.ep,
                wi.eBH    = aldex_out\$wi.eBH
            )
            aldex_df <- add_tax_label(aldex_df)
            write.table(aldex_df[order(aldex_df\$we.eBH), ],
                        file.path(out_dir, "aldex2_results.tsv"),
                        sep="\\t", quote=FALSE, row.names=FALSE, na="")

            # ALDEx2 effect plot
            pdf(file.path(out_dir, "aldex2_effect_plot.pdf"), width=8, height=6)
            aldex.plot(aldex_out, type="MA", test="welch",
                       main=paste("ALDEx2 Effect Plot:", marker))
            dev.off()

            # Volcano equivalent (effect vs -log p)
            aldex_df_nona <- aldex_df[!is.na(aldex_df\$we.eBH), ]
            p_aldex <- ggplot(aldex_df_nona, aes(x=effect, y=-log10(we.eBH+1e-300))) +
                geom_point(aes(color=we.eBH < 0.05), alpha=0.5, size=1.5) +
                scale_color_manual(values=c("FALSE"="grey60","TRUE"="#8E44AD")) +
                geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey40") +
                geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
                geom_text_repel(
                    data=subset(aldex_df_nona, we.eBH < 0.05 & abs(effect) > 1),
                    aes(label=tax_label), size=2.5, max.overlaps=15
                ) +
                theme_bw(base_size=12) +
                labs(title=paste("ALDEx2:", grp_levels[2], "vs", grp_levels[1], "-", marker),
                     x="Effect Size", y="-log10(BH p-value)")
            ggsave(file.path(out_dir, "aldex2_volcano.pdf"), p_aldex, width=9, height=7)
            ggsave(file.path(out_dir, "aldex2_volcano.png"), p_aldex, width=9, height=7, dpi=150)
        }, error=function(e) message("ALDEx2 failed: ", conditionMessage(e)))
    }

    # ── Heatmap of top differentially abundant taxa ───────────────────────
    if (length(deseq_all) > 0) {
        top_asvs <- unique(unlist(lapply(deseq_all, function(d) {
            sig <- d[!is.na(d\$padj) & d\$padj < 0.05 & abs(d\$log2FC) > 1, ]
            head(sig\$asv_id[order(sig\$padj)], 20)
        })))
        if (length(top_asvs) >= 2) {
            norm_counts <- counts(dds, normalized=TRUE)[top_asvs, , drop=FALSE]
            scaled      <- t(scale(t(norm_counts + 1)))

            tax_labels <- if ("Genus" %in% colnames(tax_tab)) {
                paste0(tax_tab[top_asvs, "Genus"], " (", top_asvs, ")")
            } else {
                top_asvs
            }
            rownames(scaled) <- tax_labels

            col_ann <- meta[colnames(scaled), grp, drop=FALSE]

            pdf(file.path(out_dir, "deseq2_significant_heatmap.pdf"), width=12, height=max(6, length(top_asvs)*0.35))
            heatmap(as.matrix(scaled),
                    margins=c(10,20),
                    col=colorRampPalette(c("#2C7BB6","white","#D7191C"))(100),
                    main=paste(marker, "- Differential Taxa (DESeq2 adj.p < 0.05)"),
                    cexRow=0.7, cexCol=0.8)
            dev.off()
        }
    }

    message("Differential abundance analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    DESeq2: ',  packageVersion('DESeq2')),
        paste0('    ALDEx2: ',  packageVersion('ALDEx2')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
