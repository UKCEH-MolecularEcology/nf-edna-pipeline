#!/usr/bin/env python3
"""
Collapse the CLARE-database ASV-level abundance table to the lowest
resolved taxon label per ASV, summing abundance for ASVs that share a
label. Ported from 12S-edna-dada2-tapirs-workflow's
05_make_cleanup_input_from_CLARE.R.

The output's "species" column is really "lowest available taxon label"
(species if resolved, else genus, else family, ... else the raw ASV ID) --
kept as "species" to match the naming the downstream blank-cleanup step
expects.

Usage:
    sintax_clare_abundance.py --in asv_taxonomy_abundance_CLARE.csv --out ncl_matrix_raw.csv
"""

import argparse
import sys

import pandas as pd

RANK_PRIORITY = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
META_COLS = ["ASV", "taxonomy"] + [c for r in RANK_PRIORITY for c in (r, f"{r}_conf")]


def taxon_label(row):
    for rank in RANK_PRIORITY:
        val = row.get(rank)
        if val and str(val).strip():
            return val
    return row["ASV"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", dest="out_path", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.in_path, dtype=str)
    for rank in RANK_PRIORITY:
        if rank in df.columns:
            df[rank] = df[rank].replace("", pd.NA)

    sample_cols = [c for c in df.columns if c not in META_COLS]
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    df["taxon_label"] = df.apply(taxon_label, axis=1)

    agg = {c: "sum" for c in sample_cols}
    agg["taxonomy"] = lambda s: next((v for v in s if pd.notna(v) and str(v).strip()), "")

    collapsed = df.groupby("taxon_label", as_index=False).agg(agg)
    collapsed = collapsed.rename(columns={"taxon_label": "species"})
    collapsed = collapsed[["species", "taxonomy"] + sample_cols]

    for c in sample_cols:
        collapsed[c] = collapsed[c].astype(int)

    collapsed.to_csv(args.out_path, index=False)


if __name__ == "__main__":
    sys.exit(main())
