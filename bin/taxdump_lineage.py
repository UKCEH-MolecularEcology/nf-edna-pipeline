#!/usr/bin/env python3
"""
Resolve NCBI taxids to a fixed 7-rank taxonomy string using a local taxdump
(nodes.dmp / names.dmp / merged.dmp), and annotate BLAST or Kraken2 output
with it. Shared by the tapirs_blast_lca and tapirs_kraken2_taxonomy modules
so the taxdump-parsing logic exists in exactly one place.

--mode coidb is a different, taxdump-free path for self-describing
references (e.g. coidb) whose BLAST hit title (stitle) already IS a
";"-joined "Kingdom;Phylum;Class;Order;Family;Genus;Species;" lineage
string -- it's parsed directly instead of resolved via an NCBI taxid, and
needs no --taxdump directory at all.

Usage:
    taxdump_lineage.py --taxdump DIR --mode blast    --in IN.blast.tsv --out OUT.blast.tax.tsv
    taxdump_lineage.py --taxdump DIR --mode kraken2  --in IN.krk       --out OUT.krk.tax.tsv
    taxdump_lineage.py             --mode coidb     --in IN.blast.tsv --out OUT.blast.tax.tsv
"""

import argparse
import os
import sys

RANKS = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")
# Presentation label used for "superkingdom" in kraken2-mode output headers/tax_rank values.
RANK_LABELS = ("domain", "phylum", "class", "order", "family", "genus", "species")
UNKNOWN = "unknown"


def import_nodes(taxdump_dir):
    taxid2parent = {}
    taxid2rank = {}
    with open(os.path.join(taxdump_dir, "nodes.dmp")) as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            if len(fields) < 3:
                continue
            try:
                taxid = int(fields[0])
                parent = int(fields[1])
            except ValueError:
                continue
            taxid2parent[taxid] = parent
            taxid2rank[taxid] = fields[2]
    return taxid2parent, taxid2rank


def import_names(taxdump_dir):
    taxid2name = {}
    with open(os.path.join(taxdump_dir, "names.dmp")) as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            if len(fields) < 4 or fields[3] != "scientific name":
                continue
            try:
                taxid = int(fields[0])
            except ValueError:
                continue
            taxid2name[taxid] = fields[1]
    return taxid2name


def import_merged(taxdump_dir):
    merged = {}
    path = os.path.join(taxdump_dir, "merged.dmp")
    if not os.path.exists(path):
        return merged
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            if len(fields) < 2:
                continue
            merged[fields[0]] = fields[1]
    return merged


class Taxdump:
    def __init__(self, taxdump_dir):
        self.taxid2parent = {}
        self.taxid2rank = {}
        self.taxid2name = {}
        self.merged = {}
        if taxdump_dir:
            self.taxid2parent, self.taxid2rank = import_nodes(taxdump_dir)
            self.taxid2name = import_names(taxdump_dir)
            self.merged = import_merged(taxdump_dir)

    def resolve_taxid(self, raw_taxid):
        """Return an int taxid, following merged.dmp if the direct lookup fails."""
        try:
            taxid = int(raw_taxid)
            if taxid in self.taxid2parent:
                return taxid
        except (ValueError, TypeError):
            pass
        new_id = self.merged.get(str(raw_taxid))
        if new_id is not None:
            try:
                taxid = int(new_id)
                if taxid in self.taxid2parent:
                    return taxid
            except ValueError:
                pass
        return None

    def get_lineage(self, taxid):
        lineage = [taxid]
        seen = {taxid}
        current = taxid
        while True:
            parent = self.taxid2parent.get(current)
            if parent is None or parent == current or parent in seen:
                break
            lineage.append(parent)
            seen.add(parent)
            current = parent
        return lineage

    def rank_name_dict(self, taxid):
        out = {}
        for node in self.get_lineage(taxid):
            rank = self.taxid2rank.get(node)
            name = self.taxid2name.get(node)
            if rank and rank != "no rank" and name:
                out[rank] = name
        return out

    def taxonomy_string(self, raw_taxid):
        """Return the 7-rank list (superkingdom..species), 'unknown' where unresolved."""
        taxid = self.resolve_taxid(raw_taxid)
        if taxid is None:
            return [UNKNOWN] * len(RANKS)

        rank_dict = self.rank_name_dict(taxid)
        values = [rank_dict.get(r, UNKNOWN) for r in RANKS]

        species = values[-1]
        genus = values[-2]
        if species != UNKNOWN:
            species = species.split("/")[0]
            tokens = species.split()
            species = "_".join(tokens[:2]) if len(tokens) >= 2 else species
            if "_sp." in species or genus == UNKNOWN:
                species = UNKNOWN
        values[-1] = species
        return values


def strip_unknown(values):
    """Drop trailing UNKNOWN entries (working backward from the deepest rank)."""
    out = list(values)
    while out and out[-1] == UNKNOWN:
        out.pop()
    return out


def run_blast_mode(taxdump, in_path, out_path):
    """Append a /-joined 7-rank taxonomy string as a 9th column to each BLAST line."""
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                fout.write(line + "\t" + "/".join([UNKNOWN] * len(RANKS)) + "\n")
                continue
            staxids = fields[3]
            if staxids in ("", "N/A") or ";" in staxids:
                tax_values = [UNKNOWN] * len(RANKS)
            else:
                tax_values = taxdump.taxonomy_string(staxids)
            fout.write(line + "\t" + "/".join(tax_values) + "\n")


def run_coidb_mode(in_path, out_path):
    """Append a /-joined 7-rank taxonomy string as a 9th column, parsed
    directly from the BLAST hit's stitle. coidb's own reference headers
    already ARE a "Kingdom;Phylum;Class;Order;Family;Genus;Species;"
    lineage, so no taxid/taxdump resolution is needed or possible (coidb
    has no real NCBI taxids)."""
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                fout.write(line + "\t" + "/".join([UNKNOWN] * len(RANKS)) + "\n")
                continue
            stitle = fields[1].strip().rstrip(";")
            parts = stitle.split(";") if stitle else []
            tax_values = [(p if p else UNKNOWN) for p in parts[:len(RANKS)]]
            tax_values += [UNKNOWN] * (len(RANKS) - len(tax_values))
            fout.write(line + "\t" + "/".join(tax_values) + "\n")


def run_kraken2_mode(taxdump, in_path, out_path):
    """Resolve each classified kraken2 read to a tax_rank/otu_id/7-rank row."""
    header = "query\ttax_rank\totu_id\t" + "\t".join(RANK_LABELS) + "\n"
    with open(out_path, "w") as fout:
        fout.write(header)
        if not os.path.exists(in_path) or os.path.getsize(in_path) == 0:
            return
        with open(in_path) as fin:
            for line in fin:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3 or fields[0] != "C":
                    continue
                seq_id, raw_taxid = fields[1], fields[2]
                if raw_taxid == "1":
                    continue
                tax_values = taxdump.taxonomy_string(raw_taxid)
                tax_str = strip_unknown(tax_values)
                if not tax_str:
                    continue
                tax_rank = RANK_LABELS[len(tax_str) - 1]
                otu_id = tax_str[-1]
                padded = tax_str + ["unidentified"] * (len(RANKS) - len(tax_str))
                fout.write(seq_id + "\t" + tax_rank + "\t" + otu_id + "\t"
                           + "\t".join(padded) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxdump", required=False, default=None,
                         help="Directory with nodes.dmp/names.dmp/merged.dmp (not used/needed for --mode coidb)")
    parser.add_argument("--mode", required=True, choices=["blast", "kraken2", "coidb"])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", dest="out_path", required=True)
    args = parser.parse_args()

    if args.mode != "coidb" and not args.taxdump:
        parser.error("--taxdump is required for --mode blast/kraken2")

    if args.mode == "coidb":
        if not os.path.exists(args.in_path) or os.path.getsize(args.in_path) == 0:
            open(args.out_path, "w").close()
            return
        run_coidb_mode(args.in_path, args.out_path)
        return

    taxdump = Taxdump(args.taxdump)

    if not os.path.exists(args.in_path) or os.path.getsize(args.in_path) == 0:
        if args.mode == "blast":
            open(args.out_path, "w").close()
        else:
            run_kraken2_mode(taxdump, args.in_path, args.out_path)
        return

    if args.mode == "blast":
        run_blast_mode(taxdump, args.in_path, args.out_path)
    else:
        run_kraken2_mode(taxdump, args.in_path, args.out_path)


if __name__ == "__main__":
    sys.exit(main())
