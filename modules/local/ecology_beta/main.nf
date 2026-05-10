process ECOLOGY_BETA {
    tag "beta_${marker}"
    label 'process_low'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/ecology:r4.3.3'

    publishDir "${params.outdir}/full_ecology/${marker}/02_beta_diversity", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata
    val group_var

    output:
    tuple val(marker), path("${marker}.beta_results/"), emit: results
    path 'versions.yml',                                 emit: versions

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
    if (length(readLines("${asv_table}")) <= 1L) {
        dir.create("${marker}.beta_results", showWarnings=FALSE, recursive=TRUE)
        writeLines("skipped: empty ASV table", "${marker}.beta_results/skipped.txt")
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


    required <- c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr")
    invisible(lapply(required, .install_pkg))
    suppressPackageStartupMessages({
        library(phyloseq); library(vegan); library(ggplot2)
        library(dplyr);    library(tidyr)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    group_var <- ${group_arg}
    out_dir   <- paste0(marker, ".beta_results")
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    set.seed(42)

    # ── Load data ────────────────────────────────────────────────────────
    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    tax_tab <- read.table("${taxonomy}",  sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    common  <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab <- asv_tab[common, , drop=FALSE]

    if (nrow(asv_tab) == 0 || ncol(asv_tab) == 0) {
        message("Empty ASV table for ", marker, " — skipping beta diversity.")
        writeLines(c('"${task.process}":', '    skipped: empty ASV table'), "versions.yml")
        writeLines("skipped: empty ASV table", file.path(out_dir, "skipped.txt"))
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

    # Relative-abundance normalisation (no rarefaction)
    ps_norm  <- transform_sample_counts(ps, function(x) x / sum(x))
    otu_mat  <- as(otu_table(ps_norm), "matrix")
    if (!taxa_are_rows(ps_norm)) otu_mat <- t(otu_mat)
    otu_t    <- t(otu_mat)   # samples x ASVs

    # CLR transformation (Aitchison distance)
    clr_transform <- function(mat) {
        mat_p <- mat + 0.5   # pseudo-count
        log(mat_p) - rowMeans(log(mat_p))
    }
    otu_clr <- clr_transform(otu_t)

    # ── 1. Distance matrices ─────────────────────────────────────────────
    dist_methods <- list(
        bray_curtis  = vegdist(otu_t,   method="bray"),
        jaccard      = vegdist(otu_t,   method="jaccard", binary=TRUE),
        aitchison    = dist(otu_clr,    method="euclidean"),
        morisita     = vegdist(otu_t,   method="morisita"),
        kulczynski   = vegdist(otu_t,   method="kulczynski"),
        horn         = vegdist(otu_t,   method="horn")
    )

    for (nm in names(dist_methods)) {
        d_mat <- as.matrix(dist_methods[[nm]])
        write.table(d_mat, file.path(out_dir, paste0(nm, "_distance.tsv")),
                    sep="\\t", quote=FALSE)
    }

    # Distance heatmap (Bray-Curtis)
    bc_mat <- as.matrix(dist_methods\$bray_curtis)
    pdf(file.path(out_dir, "bray_curtis_heatmap.pdf"), width=max(8, ncol(bc_mat)*0.3 + 2), height=max(8, ncol(bc_mat)*0.3 + 2))
    heatmap(bc_mat, symm=TRUE, col=colorRampPalette(c("white","#2C7BB6"))(100),
            margins=c(8,8), main=paste(marker, "- Bray-Curtis Distance"))
    dev.off()

    # ── 2. PERMANOVA (adonis2) ────────────────────────────────────────────
    if (!is.null(grp)) {
        meta_df  <- data.frame(sample_data(ps_norm))
        cat_vars <- names(meta_df)[sapply(meta_df, function(x) is.factor(x) || is.character(x))]
        num_vars <- names(meta_df)[sapply(meta_df, is.numeric)]
        all_vars <- c(cat_vars, num_vars)

        perm_results <- list()
        for (d_nm in c("bray_curtis","aitchison")) {
            d_obj <- dist_methods[[d_nm]]

            # Single-factor PERMANOVA for primary group variable
            formula_str <- paste0("d_obj ~ meta_df[['", grp, "']]")
            perm1 <- adonis2(as.formula(formula_str), permutations=999, parallel=${task.cpus})
            perm_results[[paste0(d_nm, "_single")]] <- data.frame(
                distance = d_nm,
                formula  = paste0("~ ", grp),
                data.frame(perm1)
            )

            # Multi-factor PERMANOVA if multiple variables available
            if (length(all_vars) > 1) {
                vars_to_use <- all_vars[seq_len(min(4, length(all_vars)))]
                tryCatch({
                    formula_multi <- paste0("d_obj ~ ", paste(vars_to_use, collapse=" + "))
                    perm_m <- adonis2(as.formula(formula_multi),
                                      data=meta_df, permutations=999, parallel=${task.cpus})
                    perm_results[[paste0(d_nm, "_multi")]] <- data.frame(
                        distance = d_nm,
                        formula  = formula_multi,
                        data.frame(perm_m)
                    )
                }, error=function(e) message("Multi-factor PERMANOVA failed: ", conditionMessage(e)))
            }
        }
        perm_df <- do.call(rbind, perm_results)
        write.table(perm_df, file.path(out_dir, "permanova_results.tsv"),
                    sep="\\t", quote=FALSE)

        # ── 3. Pairwise PERMANOVA ──────────────────────────────────────────
        groups    <- meta_df[[grp]]
        grp_levels <- unique(groups)
        if (length(grp_levels) >= 2 && length(grp_levels) <= 10) {
            pairs <- combn(grp_levels, 2, simplify=FALSE)
            pw_results <- lapply(pairs, function(pair) {
                idx <- which(groups %in% pair)
                sub_dist <- as.dist(as.matrix(dist_methods\$bray_curtis)[idx, idx])
                sub_meta <- meta_df[idx, , drop=FALSE]
                formula_str <- paste0("sub_dist ~ sub_meta[['", grp, "']]")
                pw_perm <- adonis2(as.formula(formula_str), permutations=999, parallel=${task.cpus})
                data.frame(
                    group1  = pair[1], group2 = pair[2],
                    F_value = round(pw_perm[1, "F"], 4),
                    R2      = round(pw_perm[1, "R2"], 4),
                    p_value = pw_perm[1, "Pr(>F)"]
                )
            })
            pw_df <- do.call(rbind, pw_results)
            pw_df\$p_adjusted <- p.adjust(pw_df\$p_value, method="BH")
            write.table(pw_df, file.path(out_dir, "pairwise_permanova.tsv"),
                        sep="\\t", quote=FALSE, row.names=FALSE)
        }

        # ── 4. PERMDISP (multivariate homogeneity of dispersions) ─────────
        disp_results <- list()
        for (d_nm in c("bray_curtis","aitchison")) {
            d_obj <- dist_methods[[d_nm]]
            bd    <- betadisper(d_obj, groups, type="centroid")
            bd_perm <- permutest(bd, permutations=999)
            disp_results[[d_nm]] <- data.frame(
                distance  = d_nm,
                F_value   = round(bd_perm\$tab[1,"F"], 4),
                p_value   = bd_perm\$tab[1,"Pr(>F)"],
                signif    = ifelse(bd_perm\$tab[1,"Pr(>F)"] < 0.05, "*", "ns")
            )
        }
        disp_df <- do.call(rbind, disp_results)
        write.table(disp_df, file.path(out_dir, "betadisper_permdisp.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)

        # ── 5. ANOSIM ─────────────────────────────────────────────────────
        anosim_res <- anosim(dist_methods\$bray_curtis, groups, permutations=999)
        sink(file.path(out_dir, "anosim_results.txt"))
        cat("ANOSIM (Bray-Curtis)\\n")
        cat("R statistic:", round(anosim_res\$statistic, 4), "\\n")
        cat("p-value:    ", anosim_res\$signif, "\\n\\n")
        print(anosim_res)
        sink()

        # ── 6. MRPP ───────────────────────────────────────────────────────
        mrpp_res <- mrpp(otu_t, groups, permutations=999)
        sink(file.path(out_dir, "mrpp_results.txt"))
        cat("MRPP (Bray-Curtis community distances)\\n")
        print(mrpp_res)
        sink()

        # Summary table
        summary_stats <- data.frame(
            test      = c("PERMANOVA (Bray-Curtis)", "ANOSIM (Bray-Curtis)", "MRPP"),
            statistic = c(perm_results[[1]][1,"F"], anosim_res\$statistic, mrpp_res\$delta),
            p_value   = c(perm_results[[1]][1,"Pr..F."], anosim_res\$signif, mrpp_res\$Pvalue)
        )
        summary_stats\$signif <- ifelse(summary_stats\$p_value < 0.001, "***",
                                 ifelse(summary_stats\$p_value < 0.01,  "**",
                                 ifelse(summary_stats\$p_value < 0.05,  "*", "ns")))
        write.table(summary_stats, file.path(out_dir, "multivariate_tests_summary.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
    }

    message("Beta diversity analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    vegan: ',    packageVersion('vegan')),
        paste0('    phyloseq: ', packageVersion('phyloseq')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
