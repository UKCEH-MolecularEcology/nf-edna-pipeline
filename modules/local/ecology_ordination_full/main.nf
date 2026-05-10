process ECOLOGY_ORDINATION_FULL {
    tag "ordination_full_${marker}"
    label 'process_medium'

    container 'rocker/verse:4.3.3'

    publishDir "${params.outdir}/full_ecology/${marker}/03_ordination", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata
    val group_var

    output:
    tuple val(marker), path("${marker}.ordination_results/"), emit: results
    path 'versions.yml',                                        emit: versions

    script:
    def meta_arg  = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    def group_arg = group_var ? "\"${group_var}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    # ── Package loading ──────────────────────────────────────────────────────
    r_lib <- file.path(getwd(), ".r_libs")
    dir.create(r_lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(c(r_lib, .libPaths()))

    r_ver <- numeric_version(paste(R.version\$major, R.version\$minor, sep="."))
    bioc_ver <- if (r_ver >= "4.5") "3.22" else if (r_ver >= "4.4") "3.20" else if (r_ver >= "4.3") "3.18" else "3.16"
    options(repos = c(
        BioCsoft = paste0("https://bioconductor.org/packages/", bioc_ver, "/bioc"),
        BioCann  = paste0("https://bioconductor.org/packages/", bioc_ver, "/data/annotation"),
        CRAN     = "https://cloud.r-project.org"
    ))


    required <- c("phyloseq", "vegan", "ggplot2", "dplyr", "cowplot")
    for (pkg in required) {
        if (!requireNamespace(pkg, quietly=TRUE)) {
            install.packages(pkg)
        }
    }
    suppressPackageStartupMessages({
        library(phyloseq); library(vegan); library(ggplot2)
        library(dplyr);    library(cowplot)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    group_var <- ${group_arg}
    out_dir   <- paste0(marker, ".ordination_results")
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    set.seed(42)

    # Helper: save plot
    save_plot <- function(p, name, w=9, h=7) {
        ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width=w, height=h)
        ggsave(file.path(out_dir, paste0(name, ".png")), p, width=w, height=h, dpi=150)
    }

    # ── Load data ────────────────────────────────────────────────────────
    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    tax_tab <- read.table("${taxonomy}",  sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    common  <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab <- asv_tab[common, , drop=FALSE]

    if (nrow(asv_tab) == 0 || ncol(asv_tab) == 0) {
        message("Empty ASV table for ", marker, " — skipping full ordination.")
        writeLines(c('"${task.process}":', '    skipped: empty ASV table'), "versions.yml")
        quit(status=0)
    }

    OTU <- otu_table(as.matrix(asv_tab), taxa_are_rows=TRUE)
    TAX <- tax_table(as.matrix(tax_tab[common, , drop=FALSE]))

    has_meta <- !is.null(meta_file) && file.exists(meta_file)
    if (has_meta) {
        meta     <- read.table(meta_file, sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
        common_s <- intersect(rownames(meta), colnames(asv_tab))
        meta     <- meta[common_s, , drop=FALSE]
        ps       <- phyloseq(OTU, TAX, sample_data(meta))
    } else {
        ps <- phyloseq(OTU, TAX)
    }

    grp <- if (!is.null(group_var) && has_meta && group_var %in% colnames(sample_data(ps))) {
        group_var
    } else if (has_meta) {
        colnames(sample_data(ps))[1]
    } else {
        NULL
    }

    min_reads <- min(sample_sums(ps))
    ps_rare   <- rarefy_even_depth(ps, sample.size=min_reads, rngseed=42, replace=FALSE, verbose=FALSE)
    otu_mat   <- as(otu_table(ps_rare), "matrix")
    if (!taxa_are_rows(ps_rare)) otu_mat <- t(otu_mat)
    otu_t     <- t(otu_mat)

    # CLR for Aitchison / PCA
    clr_mat <- log(otu_t + 0.5) - rowMeans(log(otu_t + 0.5))

    meta_df <- if (has_meta) data.frame(sample_data(ps_rare)) else data.frame(sample=rownames(otu_t))

    # ── Ordination helper ────────────────────────────────────────────────
    make_ord_plot <- function(scores_df, xlab, ylab, title) {
        p <- ggplot(scores_df, aes_string(x="Dim1", y="Dim2",
                                           color = if (!is.null(grp) && grp %in% colnames(scores_df)) grp else NULL,
                                           shape = if (!is.null(grp) && grp %in% colnames(scores_df)) grp else NULL)) +
            geom_point(size=3, alpha=0.85) +
            theme_bw(base_size=12) +
            labs(title=title, x=xlab, y=ylab)
        if (!is.null(grp) && grp %in% colnames(scores_df)) {
            p <- p + stat_ellipse(type="t", linetype=2, linewidth=0.5, na.rm=TRUE)
        }
        p
    }

    # ── 1. PCoA (Bray-Curtis) ────────────────────────────────────────────
    dist_bray <- vegdist(otu_t, method="bray")
    pcoa_bray <- cmdscale(dist_bray, k=min(5, nrow(otu_t)-1), eig=TRUE)
    var_exp   <- round(pcoa_bray\$eig / sum(pcoa_bray\$eig[pcoa_bray\$eig > 0]) * 100, 1)
    pcoa_df   <- data.frame(pcoa_bray\$points[,1:2], check.names=FALSE)
    colnames(pcoa_df) <- c("Dim1","Dim2")
    pcoa_df   <- cbind(pcoa_df, meta_df[rownames(pcoa_df), , drop=FALSE])

    p_pcoa <- make_ord_plot(pcoa_df,
        paste0("PCoA1 (", var_exp[1], "%)"),
        paste0("PCoA2 (", var_exp[2], "%)"),
        paste(marker, "- PCoA (Bray-Curtis)"))
    save_plot(p_pcoa, "pcoa_bray_curtis")

    write.table(data.frame(PC=paste0("PCoA",seq_along(var_exp)), Variance_pct=var_exp),
                file.path(out_dir, "pcoa_eigenvalues.tsv"), sep="\\t", quote=FALSE, row.names=FALSE)

    # ── 2. PCoA (Aitchison / CLR) ────────────────────────────────────────
    dist_ait  <- dist(clr_mat)
    pcoa_ait  <- cmdscale(dist_ait, k=min(5, nrow(clr_mat)-1), eig=TRUE)
    var_ait   <- round(pcoa_ait\$eig / sum(pcoa_ait\$eig[pcoa_ait\$eig > 0]) * 100, 1)
    pcoa_ait_df <- data.frame(pcoa_ait\$points[,1:2])
    colnames(pcoa_ait_df) <- c("Dim1","Dim2")
    pcoa_ait_df <- cbind(pcoa_ait_df, meta_df[rownames(pcoa_ait_df), , drop=FALSE])

    p_ait <- make_ord_plot(pcoa_ait_df,
        paste0("PCoA1 (", var_ait[1], "%)"),
        paste0("PCoA2 (", var_ait[2], "%)"),
        paste(marker, "- PCoA (Aitchison / CLR)"))
    save_plot(p_ait, "pcoa_aitchison")

    # ── 3. NMDS (Bray-Curtis, k=2) ───────────────────────────────────────
    nmds_fit <- metaMDS(otu_t, distance="bray", k=2, trymax=200, trace=FALSE)
    nmds_df  <- data.frame(nmds_fit\$points)
    colnames(nmds_df) <- c("Dim1","Dim2")
    nmds_df  <- cbind(nmds_df, meta_df[rownames(nmds_df), , drop=FALSE])

    p_nmds <- make_ord_plot(nmds_df,
        paste0("NMDS1"),
        paste0("NMDS2"),
        paste(marker, "- NMDS (Bray-Curtis, stress =", round(nmds_fit\$stress, 3), ")"))
    save_plot(p_nmds, "nmds_bray_curtis")

    # Shepard plot (goodness of NMDS fit)
    pdf(file.path(out_dir, "nmds_shepard.pdf"), width=7, height=6)
    stressplot(nmds_fit, main=paste(marker, "- NMDS Shepard Diagram"))
    dev.off()

    # ── 4. PCA (CLR-transformed) ─────────────────────────────────────────
    pca_res  <- prcomp(clr_mat, scale.=FALSE, center=TRUE)
    pca_var  <- round(summary(pca_res)\$importance[2,] * 100, 1)
    pca_df   <- data.frame(pca_res\$x[,1:2])
    colnames(pca_df) <- c("Dim1","Dim2")
    pca_df   <- cbind(pca_df, meta_df[rownames(pca_df), , drop=FALSE])

    p_pca <- make_ord_plot(pca_df,
        paste0("PC1 (", pca_var[1], "%)"),
        paste0("PC2 (", pca_var[2], "%)"),
        paste(marker, "- PCA (CLR-transformed)"))

    # Add top-loading taxa arrows (biplot)
    loadings   <- data.frame(pca_res\$rotation[,1:2])
    colnames(loadings) <- c("Dim1","Dim2")
    top_loads  <- loadings[order(rowSums(loadings^2), decreasing=TRUE)[1:10],]
    if ("Genus" %in% colnames(tax_tab)) {
        top_loads\$label <- tax_tab[rownames(top_loads), "Genus"]
        top_loads\$label[is.na(top_loads\$label)] <- rownames(top_loads)[is.na(top_loads\$label)]
    } else {
        top_loads\$label <- rownames(top_loads)
    }
    scale_f <- max(abs(pca_df[,c("Dim1","Dim2")])) / max(abs(top_loads[,c("Dim1","Dim2")])) * 0.7
    top_loads[,c("Dim1","Dim2")] <- top_loads[,c("Dim1","Dim2")] * scale_f

    p_pca_bi <- p_pca +
        geom_segment(data=top_loads, aes(x=0,y=0,xend=Dim1,yend=Dim2),
                     arrow=arrow(length=unit(0.2,"cm")), color="gray40", inherit.aes=FALSE) +
        geom_text(data=top_loads, aes(x=Dim1*1.1, y=Dim2*1.1, label=label),
                  size=2.5, color="gray40", inherit.aes=FALSE)
    save_plot(p_pca_bi, "pca_clr_biplot")

    write.table(data.frame(PC=paste0("PC",seq_along(pca_var)), Variance_pct=pca_var),
                file.path(out_dir, "pca_variance.tsv"), sep="\\t", quote=FALSE, row.names=FALSE)

    # ── 5. RDA with environmental variables ──────────────────────────────
    if (has_meta) {
        num_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
        cat_vars <- names(meta_df)[sapply(meta_df, function(x) is.character(x) || is.factor(x))]
        env_vars <- c(num_vars, cat_vars)

        if (length(env_vars) > 0) {
            env_df <- meta_df[rownames(otu_t), env_vars, drop=FALSE]
            env_df <- env_df[, sapply(env_df, function(x) !all(is.na(x))), drop=FALSE]

            if (ncol(env_df) > 0) {
                # Encode factors as numeric for RDA
                for (cn in cat_vars[cat_vars %in% colnames(env_df)]) {
                    env_df[[cn]] <- as.numeric(as.factor(env_df[[cn]]))
                }
                env_df <- na.omit(env_df)
                common_rows <- intersect(rownames(otu_t), rownames(env_df))

                if (length(common_rows) >= 5) {
                    # Forward selection of environmental variables
                    rda_null  <- rda(clr_mat[common_rows, ] ~ 1, data=env_df[common_rows, , drop=FALSE])
                    rda_full  <- rda(clr_mat[common_rows, ] ~ ., data=env_df[common_rows, , drop=FALSE])

                    tryCatch({
                        rda_sel  <- ordiR2step(rda_null, rda_full,
                                               permutations=199, trace=FALSE)
                        sink(file.path(out_dir, "rda_forward_selection.txt"))
                        print(summary(rda_sel))
                        sink()
                    }, error=function(e) message("RDA forward selection failed: ", conditionMessage(e)))

                    # Full RDA with all variables
                    tryCatch({
                        rda_out  <- rda(clr_mat[common_rows, ] ~ .,
                                        data=env_df[common_rows, , drop=FALSE])
                        rda_var  <- summary(rda_out)\$concont\$importance
                        rda_df   <- data.frame(scores(rda_out, display="sites"))[,1:2]
                        colnames(rda_df) <- c("Dim1","Dim2")
                        rda_df   <- cbind(rda_df, meta_df[common_rows, , drop=FALSE])
                        rda_env  <- data.frame(scores(rda_out, display="bp"))
                        colnames(rda_env)[1:2] <- c("Dim1","Dim2")
                        rda_env\$label <- rownames(rda_env)

                        scale_env <- max(abs(rda_df[,c("Dim1","Dim2")])) /
                                     max(abs(rda_env[,c("Dim1","Dim2")])) * 0.7

                        p_rda <- make_ord_plot(rda_df,
                            paste0("RDA1 (", round(rda_var[2,1]*100,1), "%)"),
                            paste0("RDA2 (", round(rda_var[2,2]*100,1), "%)"),
                            paste(marker, "- RDA (CLR)")) +
                            geom_segment(data=rda_env,
                                aes(x=0, y=0, xend=Dim1*scale_env, yend=Dim2*scale_env),
                                arrow=arrow(length=unit(0.2,"cm")), color="firebrick",
                                inherit.aes=FALSE) +
                            geom_text(data=rda_env,
                                aes(x=Dim1*scale_env*1.1, y=Dim2*scale_env*1.1, label=label),
                                size=3, color="firebrick", inherit.aes=FALSE)
                        save_plot(p_rda, "rda_environmental")

                        sink(file.path(out_dir, "rda_anova.txt"))
                        print(anova(rda_out, by="axis", permutations=999))
                        sink()
                    }, error=function(e) message("RDA failed: ", conditionMessage(e)))
                }
            }
        }
    }

    # ── 6. db-RDA (Bray-Curtis) ───────────────────────────────────────────
    if (has_meta) {
        num_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
        if (length(num_vars) >= 1) {
            env_sub  <- na.omit(meta_df[rownames(otu_t), num_vars, drop=FALSE])
            comm_r   <- intersect(rownames(otu_t), rownames(env_sub))
            if (length(comm_r) >= 5) {
                tryCatch({
                    dbrda_out <- dbrda(vegdist(otu_t[comm_r,], method="bray") ~ .,
                                       data=env_sub[comm_r, , drop=FALSE])
                    sink(file.path(out_dir, "dbrda_summary.txt"))
                    print(summary(dbrda_out))
                    print(anova(dbrda_out, by="margin", permutations=999))
                    sink()

                    dbrda_df  <- data.frame(scores(dbrda_out, display="sites"))[,1:2]
                    colnames(dbrda_df) <- c("Dim1","Dim2")
                    dbrda_df  <- cbind(dbrda_df, meta_df[comm_r, , drop=FALSE])
                    dbrda_bp  <- data.frame(scores(dbrda_out, display="bp"))
                    if (ncol(dbrda_bp) >= 2) {
                        colnames(dbrda_bp)[1:2] <- c("Dim1","Dim2")
                        dbrda_bp\$label <- rownames(dbrda_bp)
                        p_dbrda <- make_ord_plot(dbrda_df,
                            "dbRDA1", "dbRDA2",
                            paste(marker, "- db-RDA (Bray-Curtis)")) +
                            geom_segment(data=dbrda_bp,
                                aes(x=0,y=0,xend=Dim1,yend=Dim2),
                                arrow=arrow(length=unit(0.2,"cm")), color="darkgreen",
                                inherit.aes=FALSE) +
                            geom_text(data=dbrda_bp,
                                aes(x=Dim1*1.1,y=Dim2*1.1,label=label),
                                size=3, color="darkgreen", inherit.aes=FALSE)
                        save_plot(p_dbrda, "dbrda_bray_curtis")
                    }
                }, error=function(e) message("db-RDA failed: ", conditionMessage(e)))
            }
        }
    }

    # ── 7. CCA ────────────────────────────────────────────────────────────
    if (has_meta) {
        num_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
        if (length(num_vars) >= 1) {
            env_sub  <- na.omit(meta_df[rownames(otu_t), num_vars, drop=FALSE])
            comm_r   <- intersect(rownames(otu_t), rownames(env_sub))
            if (length(comm_r) >= 5) {
                tryCatch({
                    cca_out  <- cca(otu_t[comm_r, ] ~ ., data=env_sub[comm_r, , drop=FALSE])
                    sink(file.path(out_dir, "cca_summary.txt"))
                    print(summary(cca_out))
                    print(anova(cca_out, by="axis", permutations=999))
                    sink()

                    cca_df  <- data.frame(scores(cca_out, display="sites"))[,1:2]
                    colnames(cca_df) <- c("Dim1","Dim2")
                    cca_df  <- cbind(cca_df, meta_df[comm_r, , drop=FALSE])
                    p_cca   <- make_ord_plot(cca_df, "CCA1", "CCA2",
                        paste(marker, "- CCA"))
                    save_plot(p_cca, "cca_environmental")
                }, error=function(e) message("CCA failed: ", conditionMessage(e)))
            }
        }
    }

    # ── 8. Ordination panel summary ───────────────────────────────────────
    plots_to_combine <- list(p_pcoa, p_ait, p_nmds, p_pca)
    labels_panel     <- c("A) PCoA Bray-Curtis", "B) PCoA Aitchison",
                          "C) NMDS Bray-Curtis", "D) PCA CLR")
    panel <- plot_grid(plotlist=plots_to_combine, nrow=2, ncol=2,
                       labels=labels_panel, label_size=9)
    ggsave(file.path(out_dir, "ordination_panel.pdf"), panel, width=16, height=12)
    ggsave(file.path(out_dir, "ordination_panel.png"), panel, width=16, height=12, dpi=150)

    message("Full ordination analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    vegan: ',    packageVersion('vegan')),
        paste0('    phyloseq: ', packageVersion('phyloseq')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
