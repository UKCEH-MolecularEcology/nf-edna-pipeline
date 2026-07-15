#!/usr/bin/env python3
"""
Majority-vote LCA (lowest common ancestor) over BLAST hits, ported from the
12S-edna-dada2-tapirs-workflow Snakemake pipeline's mlca.py.

Input is a taxdump_lineage.py --mode blast output: 9 tab-separated columns,
no header — qseqid, stitle, sacc, staxids, pident, qcovs, evalue, bitscore,
taxonomy (7 /-joined ranks: superkingdom/phylum/class/order/family/genus/species).

For each query sequence:
  1. Keep hits with pident >= --identity and qcovs >= --coverage.
  2. Keep only hits within --bitscore percent of that query's top bitscore
     (a "top-hit bitscore window", not a fixed threshold).
  3. If fewer than --min-hits hits survive, the query is dropped entirely
     (no output row) -- its reads become "unassigned" downstream.
  4. Otherwise, vote per rank (domain..species) across the deduplicated
     taxonomy strings of the surviving hits: a rank is accepted only if its
     most common value reaches --majority percent agreement. Voting does not
     stop at the first rank that fails -- every rank is still checked -- but
     trailing "unknown" entries are stripped from the resolved chain
     afterward, so a failed deeper rank effectively truncates the result.

Usage:
    mlca.py --in IN.blast.tax.tsv --out OUT.lca.tsv \\
        --identity 98 --coverage 90 --bitscore 2 --majority 80 --min-hits 1
"""

import argparse
import sys
from collections import Counter

RANK_LABELS = ("domain", "phylum", "class", "order", "family", "genus", "species")
HEADER = "query\ttax_rank\totu_id\t" + "\t".join(RANK_LABELS) + "\tmethod\n"
UNKNOWN = "unknown"


def strip_unknown(values):
    out = list(values)
    while out and out[-1] == UNKNOWN:
        out.pop()
    return out


def load_hits(path):
    hits = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            hits.append({
                "query": fields[0],
                "pident": float(fields[4]),
                "qcovs": float(fields[5]),
                "bitscore": float(fields[7]),
                "taxonomy": fields[8],
            })
    return hits


def vote(hits_for_query, majority_frac, min_hits):
    max_bitscore = max(h["bitscore"] for h in hits_for_query)

    def process(surviving, prop):
        top = [h for h in surviving if h["bitscore"] >= max_bitscore * prop]
        if len(top) < min_hits:
            return None

        tax_strings = sorted(set(h["taxonomy"] for h in top))
        rank_lists = [t.split("/") for t in tax_strings]

        lca = []
        for rank_idx in range(len(RANK_LABELS)):
            values = [ranks[rank_idx] if rank_idx < len(ranks) else UNKNOWN for ranks in rank_lists]
            counts = Counter(values)
            winner, count = counts.most_common(1)[0]
            if count / len(tax_strings) >= majority_frac:
                lca.append(winner)

        if len(top) >= min_hits and len(top) > 1:
            method = "lca"
        elif len(top) == 1 and min_hits == 1:
            method = "single_hit"
        else:
            method = "lca"

        return lca, method

    return process


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", dest="out_path", required=True)
    parser.add_argument("--identity", type=float, required=True)
    parser.add_argument("--coverage", type=float, required=True)
    parser.add_argument("--bitscore", type=float, required=True, help="Percent bitscore drop-off window")
    parser.add_argument("--majority", type=float, required=True, help="Percent agreement required per rank")
    parser.add_argument("--min-hits", type=int, required=True, dest="min_hits")
    args = parser.parse_args()

    prop = 1 - (args.bitscore / 100.0)
    majority_frac = args.majority / 100.0

    import os
    if not os.path.exists(args.in_path) or os.path.getsize(args.in_path) == 0:
        with open(args.out_path, "w") as fout:
            fout.write(HEADER)
        return

    hits = load_hits(args.in_path)
    hits = [h for h in hits if h["pident"] >= args.identity and h["qcovs"] >= args.coverage]

    by_query = {}
    for h in hits:
        by_query.setdefault(h["query"], []).append(h)

    with open(args.out_path, "w") as fout:
        fout.write(HEADER)
        for query in sorted(by_query):
            group = by_query[query]
            result = vote(group, majority_frac, args.min_hits)(group, prop)
            if result is None:
                continue  # too few surviving top-bitscore hits -> drop the query
            lca, method = result
            lca = strip_unknown(lca)
            if lca:
                otu_id = lca[-1]
                tax_rank = RANK_LABELS[len(lca) - 1]
                padded = lca + ["unidentified"] * (len(RANK_LABELS) - len(lca))
                fout.write(f"{query}\t{tax_rank}\t{otu_id}\t" + "\t".join(padded) + f"\t{method}\n")
            else:
                padded = ["unidentified"] * len(RANK_LABELS)
                fout.write(f"{query}\tunidentified\tunidentified\t" + "\t".join(padded) + f"\t{method}\n")


if __name__ == "__main__":
    sys.exit(main())
