process ECOLOGY_ENVFIT {
    tag "envfit_${marker}"
    label 'process_low'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/ecology:r4.3.3'

    publishDir "${params.outdir}/full_ecology/${marker}/07_env_drivers", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata
    val group_var

    output:
    tuple val(marker), path("${marker}.envfit_results/"), emit: results
    path 'versions.yml',                                   emit: versions

    when:
    metadata.name != 'NO_FILE'

    script:
    def meta_arg  = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    def group_arg = group_var ? "\"${group_var}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    r_lib <- "${params.r_lib_cache}"
    dir.create(r_lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(c(r_lib, .libPaths()))
    options(mc.cores = ${task.cpus})
    if (length(readLines("${asv_table}")) <= 1L) {
        dir.create("${marker}.envfit_results", showWarnings=FALSE, recursive=TRUE)
        writeLines("skipped: empty ASV table", "${marker}.envfit_results/skipped.txt")
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

    options(repos = c(CRAN = "https://cloud.r-project.org"))


    # ── Package loading ──────────────────────────────────────────────────────
    required <- c("vegan", "ggplot2", "dplyr", "tidyr")
    invisible(lapply(required, .install_pkg))
    suppressPackageStartupMessages({
        library(vegan); library(ggplot2); library(dplyr); library(tidyr)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    group_var <- ${group_arg}
    out_dir   <- paste0(marker, ".envfit_results")
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    set.seed(42)

    if (is.null(meta_file) || !file.exists(meta_file)) {
        message("No metadata. Skipping envfit.")
        writeLines(c('"${task.process}":', '    skipped: no metadata'), "versions.yml")
        writeLines("skipped: empty ASV table", file.path(out_dir, "skipped.txt"))
        quit(status=0)
    }

    # ── Load data ────────────────────────────────────────────────────────
    asv_tab  <- read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    asv_tab[is.na(asv_tab)] <- 0
    meta     <- read.table(meta_file, sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    common_s <- intersect(colnames(asv_tab), rownames(meta))
    asv_tab  <- asv_tab[, common_s, drop=FALSE]
    meta     <- meta[common_s, , drop=FALSE]

    # Remove zero-sum samples and zero-sum ASVs
    keep_samp <- colSums(asv_tab) > 0
    asv_tab   <- asv_tab[, keep_samp, drop=FALSE]
    meta      <- meta[colnames(asv_tab), , drop=FALSE]
    keep_asv  <- rowSums(asv_tab) > 0
    asv_tab   <- asv_tab[keep_asv, , drop=FALSE]

    if (nrow(asv_tab) == 0 || ncol(asv_tab) == 0) {
        message("Empty ASV table for ", marker, " — skipping envfit analysis.")
        writeLines(c('"${task.process}":', '    skipped: empty ASV table'), "versions.yml")
        writeLines("skipped: empty ASV table", file.path(out_dir, "skipped.txt"))
        quit(status=0)
    }

    otu_t    <- t(asv_tab)
    clr_mat  <- log(otu_t + 0.5) - rowMeans(log(otu_t + 0.5))

    num_vars <- names(meta)[sapply(meta, is.numeric)]
    cat_vars <- names(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
    env_vars <- c(num_vars, cat_vars)

    if (length(env_vars) == 0) {
        message("No environmental variables found in metadata.")
        writeLines(c('"${task.process}":', '    skipped: no env vars'), "versions.yml")
        writeLines("skipped: empty ASV table", file.path(out_dir, "skipped.txt"))
        quit(status=0)
    }

    env_df   <- meta[, env_vars, drop=FALSE]
    env_clean <- env_df
    for (cn in cat_vars[cat_vars %in% colnames(env_clean)]) {
        env_clean[[cn]] <- as.numeric(as.factor(env_clean[[cn]]))
    }
    env_clean <- na.omit(env_clean)
    common_r  <- intersect(rownames(otu_t), rownames(env_clean))
    if (length(common_r) < 4) {
        message("Too few samples with complete env data (", length(common_r), ").")
        writeLines(c('"${task.process}":', '    skipped: too few samples'), "versions.yml")
        writeLines("skipped: empty ASV table", file.path(out_dir, "skipped.txt"))
        quit(status=0)
    }

    otu_sub <- otu_t[common_r, , drop=FALSE]
    env_sub <- env_clean[common_r, , drop=FALSE]

    # ── 1. envfit: vector fitting onto NMDS ordination ────────────────────
    nmds_fit <- tryCatch(
        metaMDS(otu_sub, distance="bray", k=2, trymax=200, trace=FALSE, parallel=${task.cpus}),
        error=function(e) { message("metaMDS failed: ", conditionMessage(e)); NULL }
    )

    if (is.null(nmds_fit)) {
        writeLines("skipped: metaMDS failed", file.path(out_dir, "skipped.txt"))
        writeLines(c(paste0('"${task.process}":'), '    skipped: metaMDS failed'), "versions.yml")
        quit(status=0)
    }

    ef <- tryCatch(
        envfit(nmds_fit, env_sub, permutations=999, na.rm=TRUE),
        error=function(e) { message("envfit failed: ", conditionMessage(e)); NULL }
    )

    if (!is.null(ef)) {
    sink(file.path(out_dir, "envfit_summary.txt"))
    print(ef)
    sink()
    }

    # Format envfit results
    if (!is.null(ef) && !is.null(ef\$vectors)) {
        vec_df <- data.frame(ef\$vectors\$arrows)
        vec_df\$r2      <- ef\$vectors\$r
        vec_df\$p_value <- ef\$vectors\$pvals
        vec_df\$variable <- rownames(vec_df)
        vec_df\$signif   <- ifelse(vec_df\$p_value < 0.001, "***",
                             ifelse(vec_df\$p_value < 0.01,  "**",
                             ifelse(vec_df\$p_value < 0.05,  "*", "ns")))
        vec_df <- vec_df[order(vec_df\$p_value), ]
        write.table(vec_df, file.path(out_dir, "envfit_vectors.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
    }

    # NMDS + envfit vectors plot
    nmds_df <- data.frame(nmds_fit\$points)
    colnames(nmds_df) <- c("NMDS1","NMDS2")
    nmds_df\$sample  <- rownames(nmds_df)
    grp <- if (!is.null(group_var) && group_var %in% colnames(meta[common_r,])) group_var else NULL
    if (!is.null(grp)) nmds_df[[grp]] <- meta[common_r, grp]

    # Significant env vectors
    sig_vecs <- if (!is.null(ef) && !is.null(ef\$vectors)) {
        vd <- data.frame(ef\$vectors\$arrows * sqrt(ef\$vectors\$r))
        colnames(vd) <- c("NMDS1","NMDS2")
        vd\$label    <- rownames(vd)
        vd\$p_value  <- ef\$vectors\$pvals
        vd[vd\$p_value < 0.05, , drop=FALSE]
    } else { data.frame() }

    p_ef <- ggplot(nmds_df, aes(x=NMDS1, y=NMDS2,
                                        color = if (!is.null(grp)) grp else NULL)) +
        geom_point(size=3, alpha=0.85) +
        theme_bw(base_size=12) +
        labs(title=paste(marker, "- NMDS + envfit"),
             subtitle=paste("Stress =", round(nmds_fit\$stress, 3)))

    if (nrow(sig_vecs) > 0) {
        scale_f <- max(abs(nmds_df[,c("NMDS1","NMDS2")])) /
                   max(abs(sig_vecs[,c("NMDS1","NMDS2")])) * 0.7
        p_ef <- p_ef +
            geom_segment(data=sig_vecs,
                aes(x=0, y=0, xend=NMDS1*scale_f, yend=NMDS2*scale_f),
                arrow=arrow(length=unit(0.2,"cm")), color="firebrick",
                inherit.aes=FALSE) +
            geom_text(data=sig_vecs,
                aes(x=NMDS1*scale_f*1.15, y=NMDS2*scale_f*1.15, label=label),
                size=3, color="firebrick", fontface="bold", inherit.aes=FALSE)
    }
    ggsave(file.path(out_dir, "nmds_envfit_vectors.pdf"), p_ef, width=9, height=7)
    ggsave(file.path(out_dir, "nmds_envfit_vectors.png"), p_ef, width=9, height=7, dpi=150)

    # ── 2. Mantel test ────────────────────────────────────────────────────
    num_env  <- env_sub[, num_vars[num_vars %in% colnames(env_sub)], drop=FALSE]
    if (ncol(num_env) >= 1) {
        comm_dist <- vegdist(otu_sub, method="bray")
        mantel_results <- list()
        for (v in colnames(num_env)) {
            env_vals <- num_env[[v]]
            env_dist <- dist(env_vals)
            mt <- mantel(comm_dist, env_dist, method="spearman", permutations=999, na.rm=TRUE)
            mantel_results[[v]] <- data.frame(
                variable  = v,
                mantel_r  = round(mt\$statistic, 4),
                p_value   = mt\$signif,
                signif    = ifelse(mt\$signif < 0.001, "***",
                             ifelse(mt\$signif < 0.01,  "**",
                             ifelse(mt\$signif < 0.05,  "*", "ns")))
            )
        }
        mantel_df <- do.call(rbind, mantel_results)
        write.table(mantel_df, file.path(out_dir, "mantel_tests.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)

        # Mantel plot
        p_mantel <- ggplot(mantel_df, aes(x=reorder(variable, mantel_r), y=mantel_r,
                                           fill=signif != "ns")) +
            geom_col() +
            coord_flip() +
            scale_fill_manual(values=c("FALSE"="grey70","TRUE"="#E67E22")) +
            theme_bw(base_size=12) +
            labs(title=paste(marker, "- Mantel Test (Bray-Curtis vs env)"),
                 x=NULL, y="Mantel r (Spearman)", fill="p < 0.05")
        ggsave(file.path(out_dir, "mantel_results.pdf"), p_mantel, width=8, height=max(4, ncol(num_env)*0.5))
    }

    # ── 3. Variance partitioning ──────────────────────────────────────────
    if (ncol(num_env) >= 2) {
        # Split env variables into groups for partitioning
        n_vars <- ncol(num_env)
        group1 <- num_env[, 1:ceiling(n_vars/2), drop=FALSE]
        group2 <- num_env[, (ceiling(n_vars/2)+1):n_vars, drop=FALSE]
        g1_name <- paste(colnames(group1), collapse="+")
        g2_name <- paste(colnames(group2), collapse="+")

        tryCatch({
            vp <- varpart(vegdist(otu_sub, method="bray"), group1, group2)
            sink(file.path(out_dir, "variance_partitioning.txt"))
            cat("Variance Partitioning\\n")
            cat("Group 1:", g1_name, "\\n")
            cat("Group 2:", g2_name, "\\n\\n")
            print(vp)
            sink()

            pdf(file.path(out_dir, "variance_partitioning.pdf"), width=8, height=6)
            plot(vp, digits=2, bg=c("#3498DB","#E74C3C"),
                 main=paste(marker, "- Variance Partitioning"))
            dev.off()

            vp_df <- data.frame(
                component        = c("Group1_alone","Group2_alone","Shared","Unexplained"),
                variable_group   = c(g1_name, g2_name, "Shared", "Residual"),
                r_squared        = round(c(vp\$part\$indfract[1,3],
                                           vp\$part\$indfract[2,3],
                                           vp\$part\$indfract[3,3],
                                           vp\$part\$indfract[4,3]), 4)
            )
            write.table(vp_df, file.path(out_dir, "variance_partitioning_table.tsv"),
                        sep="\\t", quote=FALSE, row.names=FALSE)
        }, error=function(e) message("Variance partitioning failed: ", conditionMessage(e)))
    } else if (ncol(num_env) >= 1) {
        message("Need >= 2 numeric env variables for partitioning. Skipping varpart.")
    }

    # ── 4. Procrustes analysis (community vs env PCA) ─────────────────────
    if (ncol(num_env) >= 2) {
        tryCatch({
            comm_pcoa <- cmdscale(vegdist(otu_sub, method="bray"),
                                  k=min(5, nrow(otu_sub)-1))
            env_pca   <- prcomp(na.omit(num_env), scale.=TRUE, center=TRUE)
            env_scores <- env_pca\$x[rownames(comm_pcoa), 1:2, drop=FALSE]

            proc_res <- protest(comm_pcoa, env_scores, permutations=999)
            sink(file.path(out_dir, "procrustes_analysis.txt"))
            cat("Procrustes Analysis: Community vs Environmental PCA\\n\\n")
            print(proc_res)
            sink()

            proc_df <- data.frame(
                metric = c("M2_statistic","correlation","p_value"),
                value  = c(round(proc_res\$ss, 4),
                           round(sqrt(1 - proc_res\$ss), 4),
                           proc_res\$signif)
            )
            write.table(proc_df, file.path(out_dir, "procrustes_stats.tsv"),
                        sep="\\t", quote=FALSE, row.names=FALSE)

            # Procrustes plot
            proc_coord <- data.frame(
                sample   = rownames(proc_res\$X),
                comm_PC1 = proc_res\$X[,1],
                comm_PC2 = proc_res\$X[,2],
                env_PC1  = proc_res\$Yrot[,1],
                env_PC2  = proc_res\$Yrot[,2]
            )
            p_proc <- ggplot(proc_coord) +
                geom_segment(aes(x=comm_PC1, y=comm_PC2, xend=env_PC1, yend=env_PC2),
                             arrow=arrow(length=unit(0.15,"cm")), color="grey50", alpha=0.7) +
                geom_point(aes(x=comm_PC1, y=comm_PC2), color="#2980B9", size=3) +
                geom_point(aes(x=env_PC1,  y=env_PC2),  color="#E74C3C", size=3, shape=17) +
                theme_bw(base_size=12) +
                labs(title=paste(marker, "- Procrustes Analysis"),
                     subtitle=paste("M2 =", round(proc_res\$ss,3),
                                    "| p =", proc_res\$signif),
                     x="Dimension 1", y="Dimension 2",
                     caption="Blue circles = community; Red triangles = environment")
            ggsave(file.path(out_dir, "procrustes_plot.pdf"), p_proc, width=8, height=7)
            ggsave(file.path(out_dir, "procrustes_plot.png"), p_proc, width=8, height=7, dpi=150)
        }, error=function(e) message("Procrustes failed: ", conditionMessage(e)))
    }

    # ── 5. Spearman correlations: alpha diversity ~ env variables ─────────
    if (ncol(num_env) >= 1) {
        alpha_obs <- specnumber(otu_sub)
        alpha_sh  <- diversity(otu_sub, index="shannon")
        alpha_df  <- data.frame(Observed=alpha_obs, Shannon=alpha_sh)

        cor_rows <- list()
        for (m in c("Observed","Shannon")) {
            for (v in colnames(num_env)) {
                ct <- cor.test(alpha_df[[m]], num_env[[v]],
                               method="spearman", exact=FALSE)
                cor_rows[[paste(m, v, sep="_")]] <- data.frame(
                    diversity_metric = m,
                    env_variable     = v,
                    rho              = round(ct\$estimate, 4),
                    p_value          = ct\$p.value,
                    signif           = ifelse(ct\$p.value < 0.001, "***",
                                       ifelse(ct\$p.value < 0.01,  "**",
                                       ifelse(ct\$p.value < 0.05,  "*", "ns")))
                )
            }
        }
        alpha_env_cor <- do.call(rbind, cor_rows)
        write.table(alpha_env_cor, file.path(out_dir, "alpha_env_correlations.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)

        # Heatmap of correlations
        cor_wide <- reshape(alpha_env_cor[, c("diversity_metric","env_variable","rho")],
                            idvar="diversity_metric", timevar="env_variable", direction="wide")
        rownames(cor_wide) <- cor_wide\$diversity_metric
        cor_wide\$diversity_metric <- NULL
        colnames(cor_wide) <- gsub("rho[.]", "", colnames(cor_wide))

        if (ncol(cor_wide) >= 2 && nrow(cor_wide) >= 1) {
            pdf(file.path(out_dir, "alpha_env_correlation_heatmap.pdf"), width=max(6, ncol(cor_wide)+2), height=4)
            heatmap(as.matrix(cor_wide),
                    col=colorRampPalette(c("#2C7BB6","white","#D7191C"))(100),
                    margins=c(10,8),
                    main=paste(marker, "- Alpha Diversity × Environmental Correlations"),
                    Rowv=NA, cexRow=0.9, cexCol=0.8)
            dev.off()
        }
    }

    message("Environmental fitting analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    vegan: ', packageVersion('vegan')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
