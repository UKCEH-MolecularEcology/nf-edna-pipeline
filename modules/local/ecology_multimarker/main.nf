process ECOLOGY_MULTIMARKER {
    tag "multimarker"
    label 'process_low'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/ecology:r4.3.3'

    publishDir "${params.outdir}/full_ecology/cross_marker", mode: 'copy'

    input:
    path asv_tables   // multiple TSVs, one per marker — staged as flat list
    path taxonomies   // matching taxonomy TSVs
    val markers       // list of marker names matching file order
    path metadata

    output:
    path "cross_marker_results/", emit: results
    path 'versions.yml',          emit: versions

    script:
    def marker_str = markers instanceof List ? markers.join(',') : markers
    def meta_arg   = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

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

    options(repos = c(CRAN = "https://cloud.r-project.org"))


    # ── Package loading ──────────────────────────────────────────────────────
    required <- c("vegan", "ggplot2", "dplyr", "tidyr", "cowplot")
    # Cross-process mutex: see ecology_alpha/main.nf for rationale.
    .lock_dir <- file.path(r_lib, ".install.lock")
    .acquired <- FALSE
    for (.i in 1:600) {
        if (dir.create(.lock_dir, showWarnings = FALSE)) { .acquired <- TRUE; break }
        Sys.sleep(1)
    }
    if (!.acquired) stop("Could not acquire R package install lock: ", .lock_dir)
    # finally (not on.exit): see ecology_alpha/main.nf for rationale.
    tryCatch({
        invisible(lapply(required, .install_pkg))
        suppressPackageStartupMessages({
            library(vegan); library(ggplot2); library(dplyr)
            library(tidyr); library(cowplot)
        })
    }, finally = { unlink(.lock_dir, recursive = TRUE) })

    out_dir    <- "cross_marker_results"
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

    meta_file  <- ${meta_arg}
    marker_str <- "${marker_str}"
    markers    <- strsplit(marker_str, ",")[[1]]
    set.seed(42)

    # ── Load all marker ASV tables ────────────────────────────────────────
    asv_files <- list.files(".", pattern="merged_asv_table[.]tsv\$", full.names=TRUE)
    tax_files <- list.files(".", pattern="[.]taxonomy[.]tsv\$",      full.names=TRUE)

    # Match files to marker names (files are staged with original names)
    data_list <- list()
    for (i in seq_along(markers)) {
        m <- markers[i]
        af <- asv_files[grep(m, asv_files, ignore.case=TRUE)]
        tf <- tax_files[grep(m, tax_files, ignore.case=TRUE)]
        if (length(af) == 0) {
            message("No ASV table found for marker: ", m, ". Skipping.")
            next
        }
        asv_t <- read.table(af[1], sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
        if (nrow(asv_t) == 0 || ncol(asv_t) == 0) {
            message("Empty ASV table for marker: ", m, ". Skipping.")
            next
        }
        tax_t <- if (length(tf) > 0) {
            read.table(tf[1], sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
        } else {
            data.frame()
        }
        data_list[[m]] <- list(asv=asv_t, tax=tax_t)
    }

    if (length(data_list) < 2) {
        message("Need at least 2 markers for cross-marker analysis. Found: ", length(data_list))
        writeLines(c('"${task.process}":', '    skipped: <2 markers'), "versions.yml")
        writeLines("skipped: empty ASV table", file.path(out_dir, "skipped.txt"))
        quit(status=0)
    }

    has_meta <- !is.null(meta_file) && file.exists(meta_file)
    meta     <- if (has_meta) {
        m <- read.table(meta_file, sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
        m
    } else {
        NULL
    }

    # Common samples across all markers
    all_samples <- lapply(data_list, function(d) colnames(d\$asv))
    common_samp <- Reduce(intersect, all_samples)
    message("Common samples across markers: ", length(common_samp))

    # ── 1. Alpha diversity comparison across markers ───────────────────────
    alpha_list <- lapply(names(data_list), function(m) {
        asv_t <- data_list[[m]]\$asv[, common_samp, drop=FALSE]
        otu_t <- t(asv_t)
        data.frame(
            sample          = rownames(otu_t),
            marker          = m,
            richness        = specnumber(otu_t),
            shannon         = diversity(otu_t, index="shannon"),
            total_reads     = rowSums(otu_t),
            n_asvs_detected = specnumber(otu_t)
        )
    })
    alpha_df <- do.call(rbind, alpha_list)
    write.table(alpha_df, file.path(out_dir, "crossmarker_alpha_diversity.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # Richness comparison plot
    p_rich <- ggplot(alpha_df, aes(x=marker, y=richness, fill=marker)) +
        geom_boxplot(alpha=0.7, outlier.shape=NA) +
        geom_jitter(width=0.15, size=1.5) +
        theme_bw(base_size=12) +
        labs(title="ASV Richness by Marker", x=NULL, y="Observed ASV Richness") +
        theme(legend.position="none")

    p_shan <- ggplot(alpha_df, aes(x=marker, y=shannon, fill=marker)) +
        geom_boxplot(alpha=0.7, outlier.shape=NA) +
        geom_jitter(width=0.15, size=1.5) +
        theme_bw(base_size=12) +
        labs(title="Shannon Diversity by Marker", x=NULL, y="Shannon H'") +
        theme(legend.position="none")

    panel_alpha <- plot_grid(p_rich, p_shan, nrow=1, align="h")
    ggsave(file.path(out_dir, "alpha_diversity_by_marker.pdf"), panel_alpha, width=12, height=6)
    ggsave(file.path(out_dir, "alpha_diversity_by_marker.png"), panel_alpha, width=12, height=6, dpi=150)

    # ── 2. Procrustes analysis between marker pairs ───────────────────────
    marker_pairs <- combn(names(data_list), 2, simplify=FALSE)
    proc_results <- list()

    for (pair in marker_pairs) {
        m1 <- pair[1]; m2 <- pair[2]
        asv1 <- t(data_list[[m1]]\$asv[, common_samp, drop=FALSE])
        asv2 <- t(data_list[[m2]]\$asv[, common_samp, drop=FALSE])

        dist1 <- vegdist(asv1, method="bray")
        dist2 <- vegdist(asv2, method="bray")

        pcoa1 <- cmdscale(dist1, k=min(4, nrow(asv1)-1))
        pcoa2 <- cmdscale(dist2, k=min(4, nrow(asv2)-1))
        k     <- min(ncol(pcoa1), ncol(pcoa2))

        tryCatch({
            proc_res <- protest(pcoa1[,1:k], pcoa2[,1:k], permutations=999)
            proc_results[[paste0(m1,"_vs_",m2)]] <- data.frame(
                marker1     = m1,
                marker2     = m2,
                M2          = round(proc_res\$ss, 4),
                correlation = round(sqrt(1 - proc_res\$ss), 4),
                p_value     = proc_res\$signif
            )
        }, error=function(e) message("Procrustes ", m1, " vs ", m2, " failed: ", conditionMessage(e)))
    }

    if (length(proc_results) > 0) {
        proc_df <- do.call(rbind, proc_results)
        write.table(proc_df, file.path(out_dir, "procrustes_between_markers.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)

        # Heatmap of correlations between markers
        if (nrow(proc_df) > 1) {
            all_ms  <- unique(c(proc_df\$marker1, proc_df\$marker2))
            cor_mat <- matrix(NA, nrow=length(all_ms), ncol=length(all_ms),
                              dimnames=list(all_ms, all_ms))
            diag(cor_mat) <- 1
            for (i in seq_len(nrow(proc_df))) {
                cor_mat[proc_df\$marker1[i], proc_df\$marker2[i]] <- proc_df\$correlation[i]
                cor_mat[proc_df\$marker2[i], proc_df\$marker1[i]] <- proc_df\$correlation[i]
            }
            pdf(file.path(out_dir, "marker_correlation_heatmap.pdf"), width=7, height=6)
            heatmap(cor_mat, symm=TRUE,
                    col=colorRampPalette(c("white","#2C7BB6"))(100),
                    main="Procrustes Correlation Between Markers",
                    margins=c(8,8))
            dev.off()
        }
    }

    # ── 3. Shared taxa (genus level) across markers ────────────────────────
    genus_by_marker <- lapply(names(data_list), function(m) {
        tax <- data_list[[m]]\$tax
        if (!"Genus" %in% colnames(tax)) return(character(0))
        unique(na.omit(tax\$Genus))
    })
    names(genus_by_marker) <- names(data_list)
    genus_by_marker <- Filter(function(x) length(x) > 0, genus_by_marker)

    if (length(genus_by_marker) >= 2) {
        # Jaccard similarity between marker genus sets
        gen_jacc <- matrix(NA, length(genus_by_marker), length(genus_by_marker),
                           dimnames=list(names(genus_by_marker), names(genus_by_marker)))
        for (i in seq_along(genus_by_marker)) {
            for (j in seq_along(genus_by_marker)) {
                a <- genus_by_marker[[i]]; b <- genus_by_marker[[j]]
                gen_jacc[i,j] <- length(intersect(a,b)) / length(union(a,b))
            }
        }
        write.table(gen_jacc, file.path(out_dir, "genus_jaccard_between_markers.tsv"),
                    sep="\\t", quote=FALSE)

        # Shared genus summary
        all_genera  <- unique(unlist(genus_by_marker))
        shared_mat  <- sapply(genus_by_marker, function(g) all_genera %in% g)
        rownames(shared_mat) <- all_genera
        n_markers_detected   <- rowSums(shared_mat)
        shared_df <- data.frame(
            genus      = all_genera,
            n_markers  = n_markers_detected,
            shared_mat,
            check.names=FALSE
        )
        shared_df <- shared_df[order(-shared_df\$n_markers), ]
        write.table(shared_df, file.path(out_dir, "genus_sharing_across_markers.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)

        # Bar chart of how many markers detect each genus
        shared_summary <- data.frame(table(n_markers_detected))
        p_shared <- ggplot(shared_summary, aes(x=n_markers_detected, y=Freq)) +
            geom_col(fill="#8E44AD", alpha=0.85) +
            theme_bw(base_size=12) +
            labs(title="Genera Shared Across Markers",
                 x="Number of markers detecting genus", y="Count of genera")
        ggsave(file.path(out_dir, "genus_sharing_barplot.pdf"), p_shared, width=7, height=5)
    }

    # ── 4. Combined ordination (all markers, all samples) ────────────────
    # Relative abundance per marker, then combine columns
    rel_list <- lapply(names(data_list), function(m) {
        asv_t    <- data_list[[m]]\$asv[, common_samp, drop=FALSE]
        rel      <- sweep(asv_t, 2, colSums(asv_t), "/")
        rownames(rel) <- paste0(m, ":", rownames(rel))
        rel
    })
    combined_asv <- do.call(rbind, rel_list)
    combined_t   <- t(combined_asv)

    dist_comb <- vegdist(combined_t, method="bray")
    pcoa_comb <- cmdscale(dist_comb, k=min(5, nrow(combined_t)-1), eig=TRUE)
    var_exp   <- round(pcoa_comb\$eig / sum(pcoa_comb\$eig[pcoa_comb\$eig>0]) * 100, 1)
    comb_df   <- data.frame(pcoa_comb\$points[,1:2])
    colnames(comb_df) <- c("PC1","PC2")
    comb_df\$sample  <- rownames(comb_df)

    first_cat <- NULL
    if (has_meta) {
        meta_sub <- meta[rownames(comb_df), , drop=FALSE]
        first_cat <- names(meta_sub)[sapply(meta_sub, function(x) is.character(x)||is.factor(x))][1]
        if (!is.null(first_cat)) comb_df[[first_cat]] <- meta_sub[[first_cat]]
    }

    p_comb <- ggplot(comb_df, aes(x=PC1, y=PC2, label=sample)) +
        geom_point(size=3, alpha=0.85,
                   aes(color=if (has_meta && !is.null(first_cat)) .data[[first_cat]] else NULL)) +
        theme_bw(base_size=12) +
        labs(title="Combined Multi-marker PCoA (Bray-Curtis)",
             x=paste0("PC1 (", var_exp[1], "%)"),
             y=paste0("PC2 (", var_exp[2], "%)"))
    ggsave(file.path(out_dir, "combined_multimarker_pcoa.pdf"), p_comb, width=9, height=7)
    ggsave(file.path(out_dir, "combined_multimarker_pcoa.png"), p_comb, width=9, height=7, dpi=150)

    message("Cross-marker analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    vegan: ',   packageVersion('vegan')),
        paste0('    ggplot2: ', packageVersion('ggplot2')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
