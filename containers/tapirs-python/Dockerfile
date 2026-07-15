FROM python:3.11-slim

LABEL org.opencontainers.image.source="https://github.com/UKCEH-MolecularEcology/nf-edna-pipeline"
LABEL org.opencontainers.image.description="Python + pandas for the 12S Tapirs/SINTAX branch (taxdump lineage lookup, MLCA, OTU table building, SINTAX parsing, blank cleanup)"
LABEL org.opencontainers.image.version="1.0"

# The scripts themselves (bin/taxdump_lineage.py, mlca.py, etc.) are NOT
# baked in here -- Nextflow automatically adds the pipeline's own bin/
# directory to PATH for every process, so they're picked up straight from
# the repo at runtime without needing a container rebuild when they change.
RUN pip install --no-cache-dir pandas==2.2.2
