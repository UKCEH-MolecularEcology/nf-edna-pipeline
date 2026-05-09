process ECOLOGY_DIVERSITY {
    tag "diversity_${marker}"
    label 'process_medium'

    container 'ghcr.io/rocker-project/verse:4.3.3'

    publishDir "${params.outdir}/ecology/${marker}/diversity", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    path metadata

    output:
    tuple val(marker), path('*.diversity_results/'), emit: results
    path 'versions.yml',                              emit: versions

    script:
    def meta_arg = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    r_lib <- file.path(getwd(), ".r_libs")
    dir.create(r_lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(c(r_lib, .libPaths()))

    r_ver <- numeric_version(paste(R.version$major, R.version$minor, sep="."))
    bioc_ver <- if (r_ver >= "4.5") "3.22" else if (r_ver >= "4.4") "3.20" else if (r_ver >= "4.3") "3.18" else "3.16"
    options(repos = c(
        BioCsoft = paste0("https://bioconductor.org/packages/", bioc_ver, "/bioc"),
        BioCann  = paste0("https://bioconductor.org/packages/", bioc_ver, "/data/annotation"),
        CRAN     = "https://cloud.r-project.org"
    ))


    pkgs <- c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr",
              "iNEXT", "microbiome", "DESeq2")
    for (pkg in pkgs) {
        if (!requireNamespace(pkg, quietly=TRUE)) {
            install.packages(pkg)
        }
    }
    suppressPackageStartupMessages({
        library(phyloseq)
        library(vegan)
        library(ggplot2)
        library(dplyr)
        library(tidyr)
    })

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    out_dir   <- paste0(marker, ".diversity_results")
    dir.create(out_dir, showWarnings = FALSE)

    # ── Load data ────────────────────────────────────────────────────────────
    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE,
                          row.names=1, check.names=FALSE)
    tax_tab <- read.table("${taxonomy}", sep="\\t", header=TRUE,
                          row.names=1, check.names=FALSE)

    # Align ASVs between tables
    common_asvs <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab     <- asv_tab[common_asvs, , drop=FALSE]
    tax_tab     <- tax_tab[common_asvs, , drop=FALSE]

    # ── Build phyloseq object ─────────────────────────────────────────────
    OTU  <- otu_table(as.matrix(asv_tab), taxa_are_rows = TRUE)
    TAX  <- tax_table(as.matrix(tax_tab))

    if (!is.null(meta_file) && file.exists(meta_file)) {
        meta    <- read.table(meta_file, sep="\\t", header=TRUE,
                              row.names=1, check.names=FALSE)
        SAMP    <- sample_data(meta)
        ps      <- phyloseq(OTU, TAX, SAMP)
    } else {
        ps      <- phyloseq(OTU, TAX)
    }

    message("Phyloseq object: ", ntaxa(ps), " taxa x ", nsamples(ps), " samples")

    # ── Alpha diversity ───────────────────────────────────────────────────
    alpha_measures <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
    alpha_df       <- estimate_richness(ps, measures = alpha_measures)
    alpha_df\$sample <- rownames(alpha_df)

    write.table(alpha_df, file.path(out_dir, "alpha_diversity.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # Alpha diversity plot
    p_alpha <- plot_richness(ps, measures = c("Observed","Shannon","Simpson")) +
        theme_bw() +
        labs(title = paste(marker, "- Alpha Diversity"),
             x = "Sample", y = "Diversity") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_dir, "alpha_diversity.pdf"), p_alpha, width=12, height=6)
    ggsave(file.path(out_dir, "alpha_diversity.png"), p_alpha, width=12, height=6, dpi=150)

    # ── Beta diversity ────────────────────────────────────────────────────
    # Rarefy to even depth for beta diversity
    min_reads <- min(sample_sums(ps))
    if (min_reads > 0) {
        ps_rare <- rarefy_even_depth(ps, sample.size = min_reads,
                                      rngseed = 42, replace = FALSE)

        # Distance matrices
        dist_bray  <- phyloseq::distance(ps_rare, method = "bray")
        dist_jacc  <- phyloseq::distance(ps_rare, method = "jaccard", binary = TRUE)

        write.table(as.matrix(dist_bray),
                    file.path(out_dir, "bray_curtis_distance.tsv"),
                    sep="\\t", quote=FALSE)
        write.table(as.matrix(dist_jacc),
                    file.path(out_dir, "jaccard_distance.tsv"),
                    sep="\\t", quote=FALSE)

        # PERMANOVA (if metadata available)
        if (!is.null(meta_file) && file.exists(meta_file)) {
            meta_df <- data.frame(sample_data(ps_rare))
            # Test first variable in metadata
            first_var <- colnames(meta_df)[1]
            perm_result <- adonis2(
                dist_bray ~ meta_df[[first_var]],
                permutations = 999
            )
            write.table(data.frame(perm_result),
                        file.path(out_dir, "permanova_bray.tsv"),
                        sep="\\t", quote=FALSE)

            # ANOSIM
            ano_result <- anosim(dist_bray, meta_df[[first_var]], permutations=999)
            sink(file.path(out_dir, "anosim_bray.txt"))
            print(ano_result)
            sink()
        }
    }

    # ── Save phyloseq object ──────────────────────────────────────────────
    saveRDS(ps, file.path(out_dir, "phyloseq_object.rds"))

    message("Diversity analysis complete. Results in: ", out_dir)

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
