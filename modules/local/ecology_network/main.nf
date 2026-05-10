process ECOLOGY_NETWORK {
    tag "network_${marker}"
    label 'process_medium'

    container 'rocker/verse:4.3.3'

    publishDir "${params.outdir}/full_ecology/${marker}/05_co_occurrence_network", mode: 'copy'

    input:
    tuple val(marker), path(asv_table), path(taxonomy)
    val min_prevalence
    val correlation_cutoff

    output:
    tuple val(marker), path("${marker}.network_results/"), emit: results
    path 'versions.yml',                                    emit: versions

    script:
    def min_prev = min_prevalence ?: 0.3
    def cor_cut  = correlation_cutoff ?: 0.6
    """
    #!/usr/bin/env Rscript

    pkgs <- c("igraph","ggraph","ggplot2","dplyr","Hmisc","psych","vegan")
    for (pkg in pkgs) {
        if (!requireNamespace(pkg, quietly=TRUE))
            install.packages(pkg, repos="https://cloud.r-project.org")
    }
    suppressPackageStartupMessages({
        library(igraph); library(ggraph); library(ggplot2)
        library(dplyr);  library(Hmisc);  library(vegan)
    })

    marker     <- "${marker}"
    min_prev   <- as.numeric("${min_prev}")
    cor_cut    <- as.numeric("${cor_cut}")
    out_dir    <- paste0(marker, ".network_results")
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    set.seed(42)

    # ── Load and filter data ─────────────────────────────────────────────
    asv_tab <- read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    tax_tab <- read.table("${taxonomy}",  sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
    common  <- intersect(rownames(asv_tab), rownames(tax_tab))
    asv_tab <- asv_tab[common, , drop=FALSE]

    # Prevalence filter: keep ASVs present in >= min_prev fraction of samples
    n_samp  <- ncol(asv_tab)
    prev    <- rowSums(asv_tab > 0) / n_samp
    asv_tab <- asv_tab[prev >= min_prev, , drop=FALSE]
    tax_tab <- tax_tab[rownames(asv_tab), , drop=FALSE]

    if (nrow(asv_tab) < 5) {
        message("Too few ASVs after prevalence filter (", nrow(asv_tab), "). Skipping network.")
        writeLines(c('"${task.process}":', '    skipped: too few taxa'), "versions.yml")
        quit(status=0)
    }

    # CLR transformation
    otu_t   <- t(asv_tab)
    otu_clr <- log(otu_t + 0.5) - rowMeans(log(otu_t + 0.5))

    # ── Correlation matrix (Spearman on CLR) ─────────────────────────────
    message("Computing Spearman correlations on ", nrow(asv_tab), " ASVs...")
    cor_mat  <- cor(otu_clr, method="spearman")
    n_obs    <- nrow(otu_clr)
    t_stat   <- cor_mat * sqrt((n_obs - 2) / (1 - cor_mat^2))
    p_mat    <- 2 * pt(-abs(t_stat), df=n_obs - 2)
    p_adj    <- matrix(p.adjust(as.vector(p_mat), method="BH"),
                       nrow=nrow(p_mat), dimnames=dimnames(p_mat))
    diag(cor_mat) <- 0

    write.table(cor_mat, file.path(out_dir, "correlation_matrix.tsv"),
                sep="\\t", quote=FALSE)
    write.table(p_adj,   file.path(out_dir, "pvalue_adj_matrix.tsv"),
                sep="\\t", quote=FALSE)

    # ── Build network with significance + strength threshold ──────────────
    cor_sig  <- cor_mat * (p_adj < 0.05) * (abs(cor_mat) >= cor_cut)
    diag(cor_sig) <- 0

    # Edge list
    edges <- which(abs(cor_sig) > 0, arr.ind=TRUE)
    edges <- edges[edges[,1] < edges[,2], , drop=FALSE]

    if (nrow(edges) == 0) {
        message("No significant edges found at threshold ", cor_cut, ". Relaxing to 0.4.")
        cor_cut_relax <- 0.4
        cor_sig <- cor_mat * (p_adj < 0.05) * (abs(cor_mat) >= cor_cut_relax)
        diag(cor_sig) <- 0
        edges   <- which(abs(cor_sig) > 0, arr.ind=TRUE)
        edges   <- edges[edges[,1] < edges[,2], , drop=FALSE]
    }

    if (nrow(edges) == 0) {
        message("Still no edges. Network analysis skipped.")
        writeLines(c('"${task.process}":', '    skipped: no edges'), "versions.yml")
        quit(status=0)
    }

    edge_df <- data.frame(
        from   = rownames(cor_sig)[edges[,1]],
        to     = colnames(cor_sig)[edges[,2]],
        weight = cor_sig[edges],
        type   = ifelse(cor_sig[edges] > 0, "positive", "negative")
    )

    write.table(edge_df, file.path(out_dir, "network_edges.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # igraph object
    g <- graph_from_data_frame(edge_df[, c("from","to","weight","type")],
                                directed=FALSE,
                                vertices=rownames(asv_tab))

    # Node attributes (taxonomy)
    V(g)\$degree      <- degree(g)
    V(g)\$betweenness <- betweenness(g, normalized=TRUE)
    V(g)\$closeness   <- closeness(g, normalized=TRUE)
    if ("Phylum" %in% colnames(tax_tab)) {
        V(g)\$phylum <- tax_tab[V(g)\$name, "Phylum"]
    } else {
        V(g)\$phylum <- "Unknown"
    }
    if ("Genus" %in% colnames(tax_tab)) {
        V(g)\$genus <- tax_tab[V(g)\$name, "Genus"]
    }

    # ── Network statistics ────────────────────────────────────────────────
    n_pos   <- sum(edge_df\$type == "positive")
    n_neg   <- sum(edge_df\$type == "negative")
    comps   <- components(g)
    net_stats <- data.frame(
        metric = c("Nodes","Edges","Positive_edges","Negative_edges",
                   "Network_density","Average_degree","Clustering_coefficient",
                   "N_components","Largest_component_size",
                   "Average_path_length","Network_diameter"),
        value  = c(
            vcount(g), ecount(g), n_pos, n_neg,
            round(graph.density(g), 4),
            round(mean(degree(g)), 2),
            round(transitivity(g, type="global"), 4),
            comps\$no,
            max(comps\$csize),
            tryCatch(round(mean_distance(g), 3), error=function(e) NA),
            tryCatch(diameter(g), error=function(e) NA)
        )
    )
    write.table(net_stats, file.path(out_dir, "network_statistics.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # ── Modularity / community detection ─────────────────────────────────
    g_pos <- subgraph.edges(g, which(E(g)\$type == "positive"), delete.vertices=FALSE)
    tryCatch({
        community <- cluster_louvain(g_pos)
        V(g)\$community <- membership(community)
        mod <- modularity(community)
        cat("Modularity (positive network, Louvain):", round(mod, 4), "\\n")
        write.table(data.frame(modularity=mod, n_modules=length(community)),
                    file.path(out_dir, "modularity.tsv"),
                    sep="\\t", quote=FALSE, row.names=FALSE)
    }, error=function(e) message("Modularity failed: ", conditionMessage(e)))

    # ── Hub taxa (high degree + betweenness) ──────────────────────────────
    node_df <- data.frame(
        asv_id       = V(g)\$name,
        degree        = V(g)\$degree,
        betweenness   = round(V(g)\$betweenness, 4),
        closeness     = round(V(g)\$closeness, 4),
        phylum        = V(g)\$phylum,
        genus         = if (!is.null(V(g)\$genus)) V(g)\$genus else NA,
        community     = if (!is.null(V(g)\$community)) V(g)\$community else NA
    )
    node_df <- node_df[order(-node_df\$degree), ]
    write.table(node_df, file.path(out_dir, "node_metrics.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # Identify keystone/hub taxa (top 5% by degree)
    hub_threshold <- quantile(node_df\$degree, 0.95)
    hub_taxa      <- node_df[node_df\$degree >= hub_threshold, ]
    write.table(hub_taxa, file.path(out_dir, "hub_keystone_taxa.tsv"),
                sep="\\t", quote=FALSE, row.names=FALSE)

    # ── Network visualization ─────────────────────────────────────────────
    # Limit to largest connected component for clarity
    main_comp <- induced_subgraph(g, V(g)[components(g)\$membership == which.max(components(g)\$csize)])

    tryCatch({
        p_net <- ggraph(main_comp, layout="fr") +
            geom_edge_link(aes(color=type, width=abs(weight)), alpha=0.4) +
            scale_edge_color_manual(values=c("positive"="#2980B9","negative"="#C0392B")) +
            scale_edge_width(range=c(0.3, 2)) +
            geom_node_point(aes(size=degree, fill=phylum),
                            shape=21, alpha=0.85) +
            scale_size(range=c(2, 10)) +
            geom_node_text(aes(label=ifelse(degree >= hub_threshold, name, "")),
                           size=2, repel=TRUE) +
            theme_graph(base_size=11) +
            labs(title=paste(marker, "- Co-occurrence Network"),
                 subtitle=paste("Nodes:", vcount(main_comp), "| Edges:", ecount(main_comp),
                               "| Spearman |r| >=", cor_cut))
        ggsave(file.path(out_dir, "network_plot.pdf"), p_net, width=12, height=10)
        ggsave(file.path(out_dir, "network_plot.png"), p_net, width=12, height=10, dpi=150)
    }, error=function(e) message("Network plot failed: ", conditionMessage(e)))

    # Positive vs negative edge pie chart
    pie_df <- data.frame(type=c("Positive","Negative"), count=c(n_pos, n_neg))
    p_pie  <- ggplot(pie_df, aes(x="", y=count, fill=type)) +
        geom_bar(stat="identity", width=1) +
        coord_polar("y") +
        scale_fill_manual(values=c("Positive"="#2980B9","Negative"="#C0392B")) +
        theme_void() +
        labs(title=paste(marker, "- Edge types"), fill="Association")
    ggsave(file.path(out_dir, "edge_type_pie.pdf"), p_pie, width=5, height=5)

    message("Network analysis complete: ", out_dir)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    igraph: ',  packageVersion('igraph')),
        paste0('    ggraph: ',  packageVersion('ggraph')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
