FROM ubuntu:22.04

LABEL org.opencontainers.image.source="https://github.com/UKCEH-MolecularEcology/nf-edna-pipeline"
LABEL org.opencontainers.image.description="usearch12 (rcedgar/usearch12, MIT-licensed) for 12S SINTAX taxonomy assignment"
LABEL org.opencontainers.image.version="12.0-beta1"

ENV USEARCH_VERSION=12.0-beta1

RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

# Pinned release binary -- never `latest` -- matching bin/install_usearch.sh
# in the source Snakemake pipeline (12S-edna-dada2-tapirs-workflow).
RUN curl -fsSL \
      "https://github.com/rcedgar/usearch12/releases/download/v${USEARCH_VERSION}/usearch_linux_x86_${USEARCH_VERSION}" \
      -o /usr/local/bin/usearch12 \
    && chmod +x /usr/local/bin/usearch12

# The SINTAX parsing script itself (bin/sintax_parse.py) is NOT baked in
# here -- Nextflow adds the pipeline's own bin/ directory to PATH for every
# process, so it's picked up straight from the repo at runtime.

RUN usearch12 --version || true
