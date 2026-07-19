#!/usr/bin/env python3
"""
Reformat a mlca.py .lca.tsv (query, tax_rank, otu_id, domain..species,
method) into the asv_id-keyed taxonomy.tsv shape build_taxonomy_table.py
already expects (asv_id, Kingdom, Phylum, Class, Order, Family, Genus,
Species, plus the LCA method/rank as extra QC columns) -- lets the existing
ASV_TAXONOMY_TABLE join logic be reused unchanged for a BLAST+LCA-derived
taxonomy, not just the DADA2 naive-Bayes one.

Usage:
    lca_to_taxonomy_table.py --in IN.lca.tsv --out OUT.taxonomy.tsv
"""

import argparse
import csv
import os
import sys

RANK_LABELS = ("domain", "phylum", "class", "order", "family", "genus", "species")
OUT_RANK_NAMES = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OUT_HEADER = ["asv_id"] + list(OUT_RANK_NAMES) + ["lca_rank", "lca_method"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", dest="out_path", required=True)
    args = parser.parse_args()

    with open(args.out_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(OUT_HEADER)

        if not os.path.exists(args.in_path) or os.path.getsize(args.in_path) == 0:
            return

        with open(args.in_path, newline="") as fin:
            reader = csv.DictReader(fin, delimiter="\t")
            for row in reader:
                out_row = [row["query"]]
                # "unidentified" is mlca.py's padding for ranks past the
                # resolved depth -- blank it out to match the DADA2
                # taxonomy.tsv convention (blank/NA = unresolved), so the
                # two taxonomy sources look consistent side by side.
                out_row += [
                    ("" if row.get(r, "") == "unidentified" else row.get(r, ""))
                    for r in RANK_LABELS
                ]
                out_row += [row.get("tax_rank", ""), row.get("method", "")]
                writer.writerow(out_row)


if __name__ == "__main__":
    sys.exit(main())
