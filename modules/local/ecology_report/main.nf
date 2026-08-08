process ECOLOGY_REPORT {
    tag "report_${marker}"
    label 'process_low'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/ecology:r4.3.3'

    publishDir "${params.outdir}/full_ecology/${marker}/00_report", mode: 'copy'

    input:
    tuple val(marker), path(alpha_dir), path(beta_dir), path(ordination_dir),
          path(network_dir), path(indicator_dir), path(envfit_dir)
    path asv_table
    path taxonomy
    path metadata

    output:
    tuple val(marker), path("${marker}_ecological_report.html"), emit: report
    path 'versions.yml',                                          emit: versions

    script:
    def meta_arg = (metadata && metadata.name != 'NO_FILE') ? "\"${metadata}\"" : 'NULL'
    """
    #!/usr/bin/env Rscript

    r_lib <- "${params.r_lib_cache}"
    dir.create(r_lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(c(r_lib, .libPaths()))
    options(mc.cores = ${task.cpus})
    if (length(readLines("${asv_table}")) <= 1L) {
        writeLines(c("<html><body><p>Skipped: empty ASV table</p></body></html>"),
                   paste0("${marker}_ecological_report.html"))
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


    pkgs <- c("rmarkdown","knitr","ggplot2","dplyr","kableExtra","DT")
    # Cross-process mutex: see ecology_alpha/main.nf for rationale.
    .lock_dir <- file.path(r_lib, ".install.lock")
    .acquired <- FALSE
    for (.i in 1:600) {
        if (dir.create(.lock_dir, showWarnings = FALSE)) { .acquired <- TRUE; break }
        Sys.sleep(1)
    }
    if (!.acquired) stop("Could not acquire R package install lock: ", .lock_dir)
    on.exit(unlink(.lock_dir, recursive = TRUE), add = TRUE)
    invisible(lapply(pkgs, .install_pkg))
    library(rmarkdown)
    unlink(.lock_dir, recursive = TRUE)

    marker    <- "${marker}"
    meta_file <- ${meta_arg}
    report_fn <- paste0(marker, "_ecological_report.Rmd")

    rmd_content <- paste0(
'---
title: "eDNA Ecological Analysis Report: ', marker, '"
date: "`r Sys.Date()`"
output:
  html_document:
    theme: flatly
    toc: true
    toc_float: true
    toc_depth: 3
    number_sections: true
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk\\$set(echo=FALSE, warning=FALSE, message=FALSE,
                        fig.align="center", out.width="100%")
library(knitr); library(dplyr)
marker <- "', marker, '"
```

# Overview {.tabset}

## Pipeline Summary

This report summarises ecological analyses for the **', marker, '** marker from the eDNA
metabarcoding pipeline. Analyses were run in the following order:

1. **Alpha Diversity** — richness, diversity, evenness, rarefaction
2. **Beta Diversity** — community dissimilarity, PERMANOVA, PERMDISP, ANOSIM
3. **Ordination** — PCoA, NMDS, PCA (CLR), RDA, db-RDA, CCA
4. **Differential Abundance** — DESeq2, ALDEx2
5. **Co-occurrence Networks** — Spearman correlation network
6. **Indicator Species** — IndVal, SIMPER, core microbiome
7. **Environmental Drivers** — envfit, Mantel, variance partitioning, Procrustes

## Data Summary

```{r data-summary}
asv_tab <- tryCatch(
  read.table("${asv_table}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE),
  error=function(e) NULL
)
tax_tab <- tryCatch(
  read.table("${taxonomy}", sep="\\t", header=TRUE, row.names=1, check.names=FALSE),
  error=function(e) NULL
)

if (!is.null(asv_tab)) {
  kable(data.frame(
    Metric = c("Number of ASVs","Number of samples",
               "Total reads","Mean reads/sample","Min reads","Max reads"),
    Value  = c(nrow(asv_tab), ncol(asv_tab),
               sum(asv_tab), round(mean(colSums(asv_tab)),0),
               min(colSums(asv_tab)), max(colSums(asv_tab)))
  ), caption="ASV table summary") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover"))
}
```

# Alpha Diversity

```{r alpha, fig.cap="Alpha diversity metrics across samples"}
alpha_file <- list.files("${alpha_dir}", pattern="alpha_diversity_metrics[.]tsv", full.names=TRUE)
if (length(alpha_file) > 0) {
  alpha_df <- read.table(alpha_file[1], sep="\\t", header=TRUE)
  kable(head(alpha_df, 20)) %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover","condensed"))
}
```

```{r alpha-plots, out.width="90%"}
plot_files <- list.files("${alpha_dir}", pattern="[.]png", full.names=TRUE)
for (f in plot_files) knitr::include_graphics(f)
```

# Beta Diversity

```{r beta-tests}
beta_file <- list.files("${beta_dir}", pattern="multivariate_tests_summary[.]tsv", full.names=TRUE)
if (length(beta_file) > 0) {
  kable(read.table(beta_file[1], sep="\\t", header=TRUE),
        caption="Multivariate community-level tests") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover"))
}
```

```{r permanova}
perm_file <- list.files("${beta_dir}", pattern="permanova_results[.]tsv", full.names=TRUE)
if (length(perm_file) > 0) {
  kable(read.table(perm_file[1], sep="\\t", header=TRUE, fill=TRUE),
        caption="PERMANOVA results") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover"))
}
```

# Ordination

```{r ordination-plots, out.width="90%", fig.cap="Ordination analyses"}
ord_files <- list.files("${ordination_dir}", pattern="[.]png", full.names=TRUE)
for (f in ord_files) knitr::include_graphics(f)
```

# Co-occurrence Network

```{r network-stats}
net_file <- list.files("${network_dir}", pattern="network_statistics[.]tsv", full.names=TRUE)
if (length(net_file) > 0) {
  kable(read.table(net_file[1], sep="\\t", header=TRUE),
        caption="Network statistics") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover"))
}
```

```{r network-plot, out.width="90%"}
net_img <- list.files("${network_dir}", pattern="network_plot[.]png", full.names=TRUE)
if (length(net_img) > 0) knitr::include_graphics(net_img[1])
```

```{r hub-taxa}
hub_file <- list.files("${network_dir}", pattern="hub_keystone_taxa[.]tsv", full.names=TRUE)
if (length(hub_file) > 0) {
  kable(read.table(hub_file[1], sep="\\t", header=TRUE),
        caption="Hub / keystone taxa (top 5% by degree)") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover"))
}
```

# Indicator Species

```{r indval}
iv_file <- list.files("${indicator_dir}", pattern="indval_significant[.]tsv", full.names=TRUE)
if (length(iv_file) > 0) {
  df <- read.table(iv_file[1], sep="\\t", header=TRUE)
  kable(head(df, 30), caption="IndVal significant indicator species (p < 0.05)") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover","condensed"))
}
```

```{r core}
core_file <- list.files("${indicator_dir}", pattern="core_microbiome_matrix[.]tsv", full.names=TRUE)
if (length(core_file) > 0) {
  core_df <- read.table(core_file[1], sep="\\t", header=TRUE)
  core50  <- core_df[core_df\$prevalence >= 50, , drop=FALSE]
  kable(head(core50[order(-core50\$prevalence),], 30),
        caption="Core microbiome (prevalence >= 50%)") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover","condensed"))
}
```

# Environmental Drivers

```{r envfit-table}
ev_file <- list.files("${envfit_dir}", pattern="envfit_vectors[.]tsv", full.names=TRUE)
if (length(ev_file) > 0) {
  kable(read.table(ev_file[1], sep="\\t", header=TRUE),
        caption="envfit: environmental vector fitting (p < 0.05 = significant)") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover"))
}
```

```{r mantel-table}
mantel_file <- list.files("${envfit_dir}", pattern="mantel_tests[.]tsv", full.names=TRUE)
if (length(mantel_file) > 0) {
  kable(read.table(mantel_file[1], sep="\\t", header=TRUE),
        caption="Mantel tests: community distance vs environmental distance") %>%
    kableExtra::kable_styling(bootstrap_options=c("striped","hover"))
}
```

```{r envfit-plots, out.width="90%"}
ev_imgs <- list.files("${envfit_dir}", pattern="[.]png", full.names=TRUE)
for (f in ev_imgs) knitr::include_graphics(f)
```

---
*Report generated by the eDNA Metabarcoding Pipeline*
*Date: `r Sys.Date()`*
')

    writeLines(rmd_content, report_fn)
    rmarkdown::render(report_fn, output_file=paste0(marker, "_ecological_report.html"),
                      quiet=TRUE)
    file.remove(report_fn)

    writeLines(c(
        paste0('"${task.process}":'),
        paste0('    rmarkdown: ', packageVersion('rmarkdown')),
        paste0('    R: ', R.version\$major, '.', R.version\$minor)
    ), "versions.yml")
    """
}
