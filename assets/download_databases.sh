#!/usr/bin/env bash
# Download reference databases for eDNA taxonomy assignment
#
# Usage:
#   bash download_databases.sh [DB_DIR] [MARKER1 MARKER2 ...]
#
# DB_DIR  : destination directory (default: databases)
# MARKERs : space-separated subset to download (default: all)
#           e.g.  bash download_databases.sh databases/ 16S CO1 12S
#
# ITS / UNITE databases require manual download — see note below.

set -euo pipefail

DB_DIR="${1:-databases}"
shift 2>/dev/null || true    # strip DB_DIR; remaining args are marker names
NEED="${*}"                  # empty = download all

mkdir -p "${DB_DIR}"
cd "${DB_DIR}"

# Returns 0 if this marker should be downloaded
need() { [[ -z "$NEED" ]] || echo " $NEED " | grep -qiw "$1"; }

# Downloads a file only if absent or empty; resumes partial downloads
fetch() {
    local url="$1" out="${2:-}"
    if [[ -n "$out" ]]; then
        [[ -s "$out" ]] && { echo "  already present: $out"; return; }
        wget -c "$url" -O "$out"
    else
        local fname
        fname="$(basename "$url")"
        [[ -s "$fname" ]] && { echo "  already present: $fname"; return; }
        wget -c "$url"
    fi
}

echo "Downloading taxonomy databases to: $(pwd)"
echo "Markers requested: ${NEED:-ALL}"
echo "========================================"

# ── 16S: SILVA (SSU, nr99, v138.1) ──────────────────────────────────────────
if need 16S; then
    echo "[16S] SILVA 138.1 DADA2-formatted"
    fetch "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
    fetch "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"
fi

# ── 18S: PR2 (v5.0.0) ───────────────────────────────────────────────────────
if need 18S; then
    echo "[18S] PR2 v5.0.0 DADA2-formatted"
    fetch "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_dada2.fasta.gz"
fi

# ── ITS: UNITE (manual) ──────────────────────────────────────────────────────
if need ITS; then
    its_file="$(ls sh_general_release_dynamic_all_*.fasta 2>/dev/null | head -1 || true)"
    if [[ -n "$its_file" ]]; then
        echo "[ITS] UNITE database already present: $its_file"
    else
        echo ""
        echo "[ITS] *** MANUAL DOWNLOAD REQUIRED ***"
        echo "  UNITE does not allow automated downloads."
        echo "  1. Go to: https://unite.ut.ee/repository.php"
        echo "  2. Download: 'UNITE general FASTA release (all eukaryotes, developer s_all)'"
        echo "  3. Place the .fasta file in: $(pwd)"
        echo "  4. Update databases.ITS.path in nextflow.config to match the filename."
        echo ""
    fi
fi

# ── CO1: MIDORI2 (GenBank 270, 2026-02-15) ──────────────────────────────────
if need CO1; then
    echo "[CO1] MIDORI2 LONGEST GB270 CO1 DADA2"
    fetch "https://www.reference-midori.info/download/Databases/GenBank270_2026-02-15/DADA2/longest/MIDORI2_LONGEST_NUC_GB270_CO1_DADA2.fasta.gz"
fi

# ── 12S: MIDORI2 (GenBank 270, 2026-02-15, srRNA) ───────────────────────────
if need 12S; then
    echo "[12S] MIDORI2 LONGEST GB270 srRNA DADA2"
    fetch "https://www.reference-midori.info/download/Databases/GenBank270_2026-02-15/DADA2/longest/MIDORI2_LONGEST_NUC_GB270_srRNA_DADA2.fasta.gz"
fi

# ── RBCL: rbcLClassifier v1 (RDP Classifier trained on NCBI plants) ─────────
if need RBCL; then
    echo "[RBCL] rbcLClassifier v1 (terrimporter)"
    fetch "https://github.com/terrimporter/rbcLClassifier/releases/download/v1/rbcLv1_trained.tar.gz"
    if [[ ! -d mydata_trained ]]; then
        echo "[RBCL] Extracting trained model..."
        tar -xzf rbcLv1_trained.tar.gz
    fi
    echo "[RBCL] Trained model: mydata_trained/"
fi

echo ""
echo "Database setup complete in: $(pwd)"
