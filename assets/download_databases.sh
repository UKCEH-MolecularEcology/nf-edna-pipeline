#!/usr/bin/env bash
# Download reference databases for eDNA taxonomy assignment
# Run once before first pipeline execution

set -euo pipefail

DB_DIR="${1:-databases}"
mkdir -p "${DB_DIR}"
cd "${DB_DIR}"

echo "Downloading taxonomy databases to: $(pwd)"
echo "========================================"

# ── 16S: SILVA (SSU, nr99, v138.1) ─────────────────────────────────────────
echo "[16S] Downloading SILVA 138.1 DADA2-formatted..."
wget -nc "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
wget -nc "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"

# ── 18S: PR2 (v5.0.0) ───────────────────────────────────────────────────────
echo "[18S] Downloading PR2 v5.0.0 DADA2-formatted..."
wget -nc "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_dada2.fasta.gz"

# ── ITS: UNITE (v10.0, all eukaryotes) ──────────────────────────────────────
echo "[ITS] Downloading UNITE v10.0..."
echo "  NOTE: UNITE requires manual download from https://unite.ut.ee/repository.php"
echo "  Download: UNITE general FASTA release (developer s_all version)"
echo "  Expected filename: sh_general_release_dynamic_s_all_*.fasta.gz"
echo "  Update databases.ITS.path in nextflow.config after downloading"

# ── CO1: MIDORI2 (GenBank 262, CO1) ─────────────────────────────────────────
echo "[CO1] Downloading MIDORI2 CO1..."
wget -nc "https://www.reference-midori.info/download/Databases/GenBank262_2024-02-02/DADA2/uniq/MIDORI2_LONGEST_NUC_SP_GB262_CO1_DADA2.fasta.gz" \
    -O MIDORI2_LONGEST_NUC_GB262_CO1_DADA2.fasta.gz || \
    echo "  Manual download from: https://www.reference-midori.info/download.e.html"

# ── 12S: MIDORI2 (GenBank 262, 12S) ─────────────────────────────────────────
echo "[12S] Downloading MIDORI2 12S..."
wget -nc "https://www.reference-midori.info/download/Databases/GenBank262_2024-02-02/DADA2/uniq/MIDORI2_LONGEST_NUC_SP_GB262_12S_DADA2.fasta.gz" \
    -O MIDORI2_LONGEST_NUC_GB262_12S_DADA2.fasta.gz || \
    echo "  Manual download from: https://www.reference-midori.info/download.e.html"

echo ""
echo "Database download complete!"
echo "Update database paths in nextflow.config if filenames differ."
ls -lh "${DB_DIR}"/*.gz 2>/dev/null || true
