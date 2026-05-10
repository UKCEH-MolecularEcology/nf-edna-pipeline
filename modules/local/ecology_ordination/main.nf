process ECOLOGY_ORDINATION {
    tag "ordination_${marker}"
    label 'process_low'

    container 'rocker/verse:4.3.3'

    publishDir "${params.outdir}/ecology/${marker}/ordination", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata

    output:
    tuple val(marker), path('*.ordination_results/'), emit: results
    path 'versions.yml',                               emit: versions

    script:
    def meta_arg = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
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


    required <- c("phyloseq", "vegan", "ggplot2", "dplyr")
    invisible(lapply(required, .install_pkg))
    suppressPackageStartupMessages({
        library(phyloseq)
        library(vegan)
        library(ggplot2)
        library(dplyr)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    out_dir   <- paste0(marker, ".ordination_results")
    dir.create(out_dir, showWarnings = FALSE)

    # Load data
    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE,
                          row.names=1, check.names=FALSE)
    tax_tab <- read.table("${taxonomy}", sep="\\t", header=TRUE,
                          row.names=1, check.names=FALSE)

    common_asvs <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab     <- asv_tab[common_asvs, , drop=FALSE]
    tax_tab     <- tax_tab[common_asvs, , drop=FALSE]

    if (nrow(asv_tab) == 0 || ncol(asv_tab) == 0) {
        message("Empty ASV table for ", marker, " — skipping ordination.")
        writeLines(c('"${task.process}":', '    skipped: empty ASV table'), "versions.yml")
        quit(status=0)
    }

    OTU <- otu_table(as.matrix(asv_tab), taxa_are_rows = TRUE)
    TAX <- tax_table(as.matrix(tax_tab))

    if (!is.null(meta_file) && file.exists(meta_file)) {
        meta  <- read.table(meta_file, sep="\\t", header=TRUE,
                            row.names=1, check.names=FALSE)
        SAMP  <- sample_data(meta)
        ps    <- phyloseq(OTU, TAX, SAMP)
        color_var <- colnames(meta)[1]
    } else {
        ps        <- phyloseq(OTU, TAX)
        color_var <- NULL
    }

    # Rarefy
    min_reads <- min(sample_sums(ps))
    if (min_reads < 1) stop("Samples with zero reads detected — check filtering.")
    ps_rare <- rarefy_even_depth(ps, sample.size=min_reads, rngseed=42, replace=FALSE)

    # ── PCoA (Bray-Curtis) ───────────────────────────────────────────────
    ord_bray <- ordinate(ps_rare, method="PCoA", distance="bray")

    p_pcoa <- plot_ordination(ps_rare, ord_bray, type="samples",
                               color = color_var) +
        geom_point(size=3) +
        stat_ellipse(aes_string(group=color_var), type="t", linetype=2, na.rm=TRUE) +
        theme_bw() +
        labs(title = paste(marker, "- PCoA (Bray-Curtis)"))
    ggsave(file.path(out_dir, "pcoa_bray_curtis.pdf"), p_pcoa, width=8, height=6)
    ggsave(file.path(out_dir, "pcoa_bray_curtis.png"), p_pcoa, width=8, height=6, dpi=150)

    # Save eigenvalues
    write.table(
        data.frame(PC=paste0("PC", seq_along(ord_bray\$values\$Eigenvalues)),
                   Eigenvalue=ord_bray\$values\$Eigenvalues,
                   Variance_explained=ord_bray\$values\$Relative_eig*100),
        file.path(out_dir, "pcoa_eigenvalues.tsv"),
        sep="\\t", quote=FALSE, row.names=FALSE
    )

    # ── NMDS (Bray-Curtis) ───────────────────────────────────────────────
    ord_nmds <- ordinate(ps_rare, method="NMDS", distance="bray",
                         trymax=100, k=2)

    p_nmds <- plot_ordination(ps_rare, ord_nmds, type="samples",
                               color = color_var) +
        geom_point(size=3) +
        annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
                 label=paste0("Stress = ", round(ord_nmds\$stress, 3))) +
        theme_bw() +
        labs(title = paste(marker, "- NMDS (Bray-Curtis)"))
    ggsave(file.path(out_dir, "nmds_bray_curtis.pdf"), p_nmds, width=8, height=6)
    ggsave(file.path(out_dir, "nmds_bray_curtis.png"), p_nmds, width=8, height=6, dpi=150)

    # ── RDA/CCA (if metadata continuous variables present) ───────────────
    if (!is.null(meta_file) && file.exists(meta_file)) {
        meta_df <- data.frame(sample_data(ps_rare))
        numeric_vars <- names(meta_df)[sapply(meta_df, is.numeric)]

        if (length(numeric_vars) > 0) {
            otu_mat <- t(otu_table(ps_rare))
            rda_out <- rda(otu_mat ~ ., data=meta_df[, numeric_vars, drop=FALSE],
                           scale=FALSE)
            sink(file.path(out_dir, "rda_summary.txt"))
            print(summary(rda_out))
            sink()

            pdf(file.path(out_dir, "rda_triplot.pdf"), width=8, height=8)
            plot(rda_out, display=c("sites","bp"), type="n",
                 main=paste(marker, "- RDA"))
            points(rda_out, display="sites", pch=16, cex=1)
            text(rda_out, display="bp", col="blue", arrow.mul=0.8)
            dev.off()
        }
    }

    message("Ordination analysis complete.")

    writeLines(
        c(
            paste0('"${task.process}":'),
            paste0('    phyloseq: ', packageVersion('phyloseq')),
            paste0('    vegan: ',    packageVersion('vegan')),
            paste0('    R: ', R.version\$major, '.', R.version\$minor)
        ),
        "versions.yml"
    )
    """
}
