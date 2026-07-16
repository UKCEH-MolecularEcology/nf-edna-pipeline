#!/usr/bin/env python3
"""
Join a marker's merged ASV abundance table with its collective taxonomy
assignment into one combined table, keyed consistently via the ASV
lookup produced by MERGE_ASV_TABLES.

Produces two files from the same join:
  --out-abundance   ASV, sequence, <sample columns...>, <taxonomy columns...>
                     (the human-facing deliverable: every ASV's sequence,
                     per-sample counts, and taxonomy lineage in one place)
  --out-by-sequence sequence, <taxonomy columns...>
                     (re-keyed by sequence instead of ASV label, no sample
                     columns -- feeds the ecology modules, which join the
                     taxonomy file to the ASV table by matching column 1
                     against merged_asv_table.tsv's sequence-keyed asv_id)

Usage:
    build_taxonomy_table.py --asv-lookup asv_lookup.tsv \\
        --merged-table merged_asv_table.tsv --taxonomy taxonomy.tsv \\
        --out-abundance MARKER.asv_taxonomy_abundance.tsv \\
        --out-by-sequence MARKER.taxonomy_by_sequence.tsv
"""

import argparse
import csv
import os
import sys


def read_tsv(path):
    """Return (header, rows) with rows as list of dicts, or ([], []) if empty/missing."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)
    return header, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--asv-lookup", required=True)
    parser.add_argument("--merged-table", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--out-abundance", required=True)
    parser.add_argument("--out-by-sequence", required=True)
    args = parser.parse_args()

    lookup_header, lookup_rows = read_tsv(args.asv_lookup)
    table_header, table_rows = read_tsv(args.merged_table)
    tax_header, tax_rows = read_tsv(args.taxonomy)

    if not lookup_rows:
        # Zero-ASV marker: write header-only outputs.
        with open(args.out_abundance, "w", newline="") as fh:
            fh.write("ASV\tsequence\n")
        with open(args.out_by_sequence, "w", newline="") as fh:
            fh.write("sequence\n")
        return

    sample_cols = [c for c in table_header if c != "asv_id"]
    counts_by_seq = {row["asv_id"]: row for row in table_rows}

    tax_cols = [c for c in tax_header if c != "asv_id"]
    tax_by_asv = {row["asv_id"]: row for row in tax_rows}

    abundance_fields = ["ASV", "sequence"] + sample_cols + tax_cols
    by_seq_fields = ["sequence"] + tax_cols

    with open(args.out_abundance, "w", newline="") as f_abund, \
         open(args.out_by_sequence, "w", newline="") as f_seq:

        w_abund = csv.DictWriter(f_abund, fieldnames=abundance_fields, delimiter="\t")
        w_abund.writeheader()
        w_seq = csv.DictWriter(f_seq, fieldnames=by_seq_fields, delimiter="\t")
        w_seq.writeheader()

        for row in lookup_rows:
            asv, seq = row["ASV"], row["sequence"]
            counts = counts_by_seq.get(seq, {})
            tax = tax_by_asv.get(asv, {})

            abund_row = {"ASV": asv, "sequence": seq}
            for c in sample_cols:
                abund_row[c] = counts.get(c, 0)
            for c in tax_cols:
                abund_row[c] = tax.get(c, "")
            w_abund.writerow(abund_row)

            seq_row = {"sequence": seq}
            for c in tax_cols:
                seq_row[c] = tax.get(c, "")
            w_seq.writerow(seq_row)


if __name__ == "__main__":
    sys.exit(main())
