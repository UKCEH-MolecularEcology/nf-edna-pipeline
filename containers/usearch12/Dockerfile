FROM ubuntu:22.04

LABEL org.opencontainers.image.source="https://github.com/UKCEH-MolecularEcology/nf-edna-pipeline"
LABEL org.opencontainers.image.description="usearch12 (rcedgar/usearch12, MIT-licensed) for 12S SINTAX taxonomy assignment"
LABEL org.opencontainers.image.version="12.0-beta1"

ENV USEARCH_RELEASE_TAG=v12.0-beta1
ENV USEARCH_ASSET_NAME=usearch_linux_x86_12.0-beta

# python3: sintax_parse.py (a bin/ script) runs inside this same container
#   right after usearch12 -sintax, and ubuntu:22.04 ships no Python by default.
# procps: provides `ps`, which Nextflow requires inside every container to
#   poll task resource usage (also missing from the bare base image).
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl python3 procps \
    && rm -rf /var/lib/apt/lists/*

# Pinned release binary -- never `latest` -- matching bin/install_usearch.sh
# in the source Snakemake pipeline (12S-edna-dada2-tapirs-workflow). Note the
# release TAG (v12.0-beta1) and the release ASSET filename (12.0-beta, no
# trailing "1") differ -- that mismatch is real, not a typo here.
RUN curl -fsSL \
      "https://github.com/rcedgar/usearch12/releases/download/${USEARCH_RELEASE_TAG}/${USEARCH_ASSET_NAME}" \
      -o /usr/local/bin/usearch12 \
    && chmod +x /usr/local/bin/usearch12

# The SINTAX parsing script itself (bin/sintax_parse.py) is NOT baked in
# here -- Nextflow adds the pipeline's own bin/ directory to PATH for every
# process, so it's picked up straight from the repo at runtime.

RUN usearch12 --version || true
