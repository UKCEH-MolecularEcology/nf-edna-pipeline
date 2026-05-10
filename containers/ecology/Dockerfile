FROM rocker/verse:4.3.3

LABEL org.opencontainers.image.source="https://github.com/UKCEH-MolecularEcology/nf-edna-pipeline"
LABEL org.opencontainers.image.description="eDNA ecology analysis: phyloseq, vegan, DESeq2, ALDEx2, microbiome"
LABEL org.opencontainers.image.version="1.0"

# R 4.3.x requires Bioconductor 3.18
ENV BIOC_VERSION=3.18

# CRAN packages (no Bioconductor dependencies)
RUN Rscript -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages( \
        c('vegan', 'cowplot', 'indicspecies', 'iNEXT', 'kableExtra', 'DT'), \
        Ncpus = parallel::detectCores() \
    )"

# Bioconductor packages (phyloseq, DESeq2, ALDEx2, microbiome)
RUN Rscript -e "\
    options(repos = c( \
        BioCsoft = paste0('https://bioconductor.org/packages/', Sys.getenv('BIOC_VERSION'), '/bioc'), \
        BioCann  = paste0('https://bioconductor.org/packages/', Sys.getenv('BIOC_VERSION'), '/data/annotation'), \
        CRAN     = 'https://cloud.r-project.org' \
    )); \
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager'); \
    BiocManager::install( \
        c('phyloseq', 'DESeq2', 'ALDEx2', 'microbiome'), \
        update = FALSE, ask = FALSE, Ncpus = parallel::detectCores() \
    )"

# Verify everything loaded — fail the build if any package is missing
RUN Rscript -e "\
    pkgs <- c( \
        'vegan', 'cowplot', 'indicspecies', 'iNEXT', 'kableExtra', 'DT', \
        'phyloseq', 'DESeq2', 'ALDEx2', 'microbiome', \
        'igraph', 'ggplot2', 'dplyr', 'tidyr', 'knitr', 'rmarkdown' \
    ); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop('Missing packages: ', paste(missing, collapse = ', ')); \
    cat('All', length(pkgs), 'packages verified.\n')"
