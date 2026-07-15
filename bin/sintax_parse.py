#!/usr/bin/env python3
"""
Parse usearch12 `-sintax -tabbedout` output into per-ASV taxonomy/confidence
columns, then build a species-level and an ASV-level abundance table.
Ported from 12S-edna-dada2-tapirs-workflow's 03a_sintax_run.R.

SINTAX tabbed output: column 1 = query ID, column 2 = a comma-separated
string of rank:name(confidence) pairs, e.g.
    k:Eukaryota(1.00),p:Chordata(0.98),...,s:Danio_rerio(0.87)
(columns 3/4, if present, are the predicted strand / bootstrap detail and
are ignored here, matching the original R parser.)

Usage:
    sintax_parse.py --sintax-out OUT.sintax.tsv --asv-lookup asv_lookup.tsv \\
        --abundance-table merged_asv_table.tsv --db-name CLARE \\
        --out-prefix marker_CLARE
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

RANK_PREFIXES = [("k", "kingdom"), ("p", "phylum"), ("c", "class"),
                  ("o", "order"), ("f", "family"), ("g", "genus"), ("s", "species")]


def extract_ranks(taxonomy_str):
    """Return {rank: name, rank_conf: value} for one SINTAX taxonomy string."""
    out = {}
    for prefix, rank in RANK_PREFIXES:
        m = re.search(rf"{prefix}:([^,(]+)\(([^)]+)\)", taxonomy_str or "")
        out[rank] = m.group(1) if m else ""
        out[f"{rank}_conf"] = m.group(2) if m else ""
    return out


def read_sintax_simple(path):
    """Read the first two columns of the tabbed SINTAX output: ASV, taxonomy."""
    rows = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return rows
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            rows[fields[0]] = fields[1]
    return rows


def read_asv_lookup(path):
    lookup = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            lookup.append(row)
    return lookup


def read_abundance_table(path):
    """merged_asv_table.tsv: first column asv_id (= ASV SEQUENCE), then one
    column per sample. Returns (samples, {sequence: {sample: count}})."""
    samples = []
    by_seq = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return samples, by_seq
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return samples, by_seq
        samples = header[1:]
        for row in reader:
            if not row:
                continue
            seq = row[0]
            by_seq[seq] = {s: int(v) for s, v in zip(samples, row[1:])}
    return samples, by_seq


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sintax-out", required=True)
    parser.add_argument("--asv-lookup", required=True)
    parser.add_argument("--abundance-table", required=True)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    sintax_rows = read_sintax_simple(args.sintax_out)
    lookup = read_asv_lookup(args.asv_lookup)
    samples, abundance_by_seq = read_abundance_table(args.abundance_table)

    parsed_fields = ["ASV", "taxonomy", "database"] + \
        [f for _, rank in RANK_PREFIXES for f in (rank, f"{rank}_conf")]

    with open(f"{args.out_prefix}_parsed.tsv", "w", newline="") as fparsed, \
         open(f"{args.out_prefix}_asv_taxonomy.tsv", "w", newline="") as ftax:

        w_parsed = csv.writer(fparsed, delimiter="\t")
        w_parsed.writerow(parsed_fields)

        tax_fields = ["ASV", "sequence", "length"] + \
            [f for _, rank in RANK_PREFIXES for f in (rank, f"{rank}_conf")]
        w_tax = csv.writer(ftax, delimiter="\t")
        w_tax.writerow(tax_fields)

        asv_species = {}       # ASV -> species (for species_abundance)
        asv_ranks = {}         # ASV -> dict of rank/conf values
        asv_to_seq = {}

        for row in lookup:
            asv = row["ASV"]
            asv_to_seq[asv] = row["sequence"]
            taxonomy_str = sintax_rows.get(asv, "")
            ranks = extract_ranks(taxonomy_str)
            asv_ranks[asv] = ranks
            if ranks["species"]:
                asv_species[asv] = ranks["species"]

            w_parsed.writerow([asv, taxonomy_str, args.db_name] +
                               [ranks[f] for _, rank in RANK_PREFIXES for f in (rank, f"{rank}_conf")])
            w_tax.writerow([asv, row["sequence"], row.get("length", "")] +
                           [ranks[f] for _, rank in RANK_PREFIXES for f in (rank, f"{rank}_conf")])

    # species_abundance.csv: samples as rows, species-collapsed ASV counts as columns
    with open(f"species_abundance_{args.db_name}.csv", "w", newline="") as fh:
        if not asv_species:
            csv.writer(fh).writerow(["sample"])
        else:
            species_list = sorted(set(asv_species.values()))
            w = csv.writer(fh)
            w.writerow(["sample"] + species_list)
            for sample in samples:
                counts = {sp: 0 for sp in species_list}
                for asv, species in asv_species.items():
                    seq = asv_to_seq.get(asv)
                    counts[species] += abundance_by_seq.get(seq, {}).get(sample, 0)
                w.writerow([sample] + [counts[sp] for sp in species_list])

    # asv_taxonomy_abundance.csv: one row per ASV, sample columns + taxonomy/conf columns
    with open(f"asv_taxonomy_abundance_{args.db_name}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        rank_cols = [f for _, rank in RANK_PREFIXES for f in (rank, f"{rank}_conf")]
        w.writerow(["ASV"] + samples + ["taxonomy"] + rank_cols)
        for row in lookup:
            asv = row["ASV"]
            seq = row["sequence"]
            counts = [abundance_by_seq.get(seq, {}).get(s, 0) for s in samples]
            ranks = asv_ranks.get(asv, {})
            taxonomy_str = sintax_rows.get(asv, "")
            w.writerow([asv] + counts + [taxonomy_str] + [ranks.get(c, "") for c in rank_cols])


if __name__ == "__main__":
    sys.exit(main())
