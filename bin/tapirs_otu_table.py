#!/usr/bin/env python3
"""
Build a sample x OTU abundance matrix (+ taxonomy string) from a collection
of per-sample `.lca.tsv` (BLAST/MLCA route) or `.krk.tax.tsv` (Kraken2 route)
files, plus each sample's rereplicated FASTA (used only as a total-read-count
denominator for the "unassigned" row). Ported from the
12S-edna-dada2-tapirs-workflow Snakemake pipeline's mlca_to_tsv.py /
mlca_to_tsv_full_lineage.py -- both flavors are emitted from one pass:

  --out-regular      #OTU_ID = the deepest resolved taxon name (or "unassigned")
  --out-full-lineage #OTU_ID = the full "d__X; p__Y; ..." lineage string itself

Each input tax file's rows are truncated to --highest-rank..--lowest-rank
(ranks above --highest-rank are dropped entirely; ranks below --lowest-rank
are collapsed away), matching the Snakemake pipeline's rank-window config.

Usage:
    tapirs_otu_table.py --tax-files S1.lca.tsv S2.lca.tsv ... \\
        --rerep-files S1.rerep.fasta S2.rerep.fasta ... \\
        --lowest-rank species --highest-rank order \\
        --out-regular OUT.tsv --out-full-lineage OUT_full_lineage.tsv
"""

import argparse
import os
import re
import sys
from collections import defaultdict

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
SIZE_RE = re.compile(r";size=(\d+)")


def sample_name(path):
    return os.path.basename(path).split(".")[0]


def count_fasta_records(path):
    if not path or not os.path.exists(path):
        return 0
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def parse_tax_file(path, min_rank_idx, max_rank_idx):
    """Return {otu_id: count} and {otu_id: rank_vals} for one sample's tax file."""
    counts = defaultdict(float)
    taxonomy = {}
    ranks_upper = RANKS[max_rank_idx:]
    lowest_rank_name = RANKS[min_rank_idx]

    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return counts, taxonomy

    with open(path) as fh:
        lines = fh.readlines()

    if len(lines) <= 1:
        return counts, taxonomy

    for line in lines[1:]:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        query, tax_rank, otu_id = fields[0], fields[1], fields[2]
        rank_vals = list(fields[3:10])

        if tax_rank not in ranks_upper:
            continue

        # Collapse anything below the configured lowest rank.
        for i in range(min_rank_idx + 1, len(RANKS)):
            rank_vals[i] = "unidentified"
        if RANKS.index(tax_rank) > min_rank_idx:
            tax_rank = lowest_rank_name

        if rank_vals[min_rank_idx] != "unidentified":
            otu_id = rank_vals[min_rank_idx]

        m = SIZE_RE.search(query)
        abundance = float(m.group(1)) if m else 1.0

        counts[otu_id] += abundance
        if otu_id not in taxonomy:
            taxonomy[otu_id] = rank_vals

    return counts, taxonomy


def lineage_string(rank_vals):
    prefixes = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")
    resolved = []
    for prefix, val in zip(prefixes, rank_vals):
        if val == "unidentified":
            break
        resolved.append(f"{prefix}{val}")
    return "; ".join(resolved)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tax-files", nargs="+", required=True)
    parser.add_argument("--rerep-files", nargs="+", required=True)
    parser.add_argument("--lowest-rank", default="species", choices=RANKS)
    parser.add_argument("--highest-rank", default="order", choices=RANKS)
    parser.add_argument("--out-regular", required=True)
    parser.add_argument("--out-full-lineage", required=True)
    args = parser.parse_args()

    min_rank_idx = RANKS.index(args.lowest_rank)
    max_rank_idx = RANKS.index(args.highest_rank)

    rerep_by_sample = {sample_name(p): p for p in args.rerep_files}

    abundance = defaultdict(dict)     # otu_id -> sample -> count
    taxonomy = {}                     # otu_id -> rank_vals
    samples = []
    unassigned = {}

    for tax_path in args.tax_files:
        sample = sample_name(tax_path)
        samples.append(sample)
        counts, tax_for_sample = parse_tax_file(tax_path, min_rank_idx, max_rank_idx)

        for otu_id, rank_vals in tax_for_sample.items():
            if otu_id not in taxonomy:
                taxonomy[otu_id] = rank_vals

        assigned_reads = sum(c for otu, c in counts.items() if otu != "unidentified")
        for otu_id, c in counts.items():
            if otu_id == "unidentified":
                continue
            abundance[otu_id][sample] = abundance[otu_id].get(sample, 0) + c

        total_reads = count_fasta_records(rerep_by_sample.get(sample))
        unassigned[sample] = max(total_reads - assigned_reads, 0)

    samples = sorted(set(samples))
    otu_ids = sorted(taxonomy.keys())

    def write_table(out_path, use_full_lineage_id):
        with open(out_path, "w") as fout:
            rows = []
            for otu_id in otu_ids:
                rank_vals = taxonomy[otu_id]
                tax_str = lineage_string(rank_vals)
                row_id = tax_str if use_full_lineage_id else otu_id
                counts = [int(abundance[otu_id].get(s, 0)) for s in samples]
                depth = tax_str.count(";") + 1 if tax_str else 0
                rows.append((depth, row_id, tax_str, counts))

            rows.sort(key=lambda r: (-r[0], r[1]))

            unassigned_id = "u__unassigned" if use_full_lineage_id else "unassigned"
            unassigned_counts = [int(unassigned.get(s, 0)) for s in samples]

            header = ["#OTU_ID"] + samples
            if not use_full_lineage_id:
                header += ["taxonomy"]
            fout.write("\t".join(header) + "\n")

            for depth, row_id, tax_str, counts in rows:
                fields = [row_id] + [str(c) for c in counts]
                if not use_full_lineage_id:
                    fields += [tax_str]
                fout.write("\t".join(fields) + "\n")

            unassigned_fields = [unassigned_id] + [str(c) for c in unassigned_counts]
            if not use_full_lineage_id:
                unassigned_fields += ["u__unassigned"]
            fout.write("\t".join(unassigned_fields) + "\n")

    write_table(args.out_regular, use_full_lineage_id=False)
    write_table(args.out_full_lineage, use_full_lineage_id=True)


if __name__ == "__main__":
    sys.exit(main())
