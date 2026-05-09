#!/usr/bin/env python3
"""
Parse RDP Classifier fixrank output to a standard taxonomy TSV.

Usage:
    parse_rdp.py <rdp_fixrank.txt> <out.taxonomy.tsv> <bootstrap_cutoff>

The cutoff (0-1) is applied to all ranks: assignments below the threshold
are left blank.  genus_boot and species_boot columns always contain the
raw confidence value regardless of the cutoff.
"""

import csv
import sys

OUT_RANKS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# Map output rank names to the RDP rank labels found in the trained model.
# The rbcLClassifier uses NCBI-sourced taxonomy; rank labels are lowercase.
RANK_LOOKUP = {
    'Kingdom': ['superkingdom', 'kingdom'],
    'Phylum':  ['phylum'],
    'Class':   ['class'],
    'Order':   ['order'],
    'Family':  ['family'],
    'Genus':   ['genus'],
    'Species': ['species'],
}


def parse_fixrank_line(cols):
    """Return (rank_val, rank_boot) dicts from a fixrank output line.

    RDP fixrank format (cols[2:]):
        rank_label  taxon_name  confidence  rank_label  taxon_name  confidence ...
    """
    rank_val = {}
    rank_boot = {}
    i = 2
    while i + 2 < len(cols):
        rank = cols[i].lower()
        taxon = cols[i + 1]
        try:
            conf = float(cols[i + 2])
        except (ValueError, IndexError):
            conf = 0.0
        rank_val[rank] = taxon
        rank_boot[rank] = conf
        i += 3
    return rank_val, rank_boot


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} rdp_raw.txt out.tsv bootstrap_cutoff")

    in_file, out_file = sys.argv[1], sys.argv[2]
    try:
        cutoff = float(sys.argv[3])
    except ValueError:
        sys.exit(f"ERROR: bootstrap_cutoff must be a number, got '{sys.argv[3]}'")

    with open(in_file) as fh, open(out_file, 'w', newline='') as out_fh:
        writer = csv.writer(out_fh, delimiter='\t')
        writer.writerow(['asv_id'] + OUT_RANKS + ['genus_boot', 'species_boot'])

        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3 or not cols[0]:
                continue

            seq_name = cols[0]
            rank_val, rank_boot = parse_fixrank_line(cols)

            row = [seq_name]
            for r in OUT_RANKS:
                val = ''
                for rk in RANK_LOOKUP[r]:
                    if rk in rank_val and rank_boot.get(rk, 0.0) >= cutoff:
                        val = rank_val[rk]
                        break
                row.append(val)

            row.append(f"{rank_boot.get('genus', ''):.4f}" if 'genus' in rank_boot else '')
            row.append(f"{rank_boot.get('species', ''):.4f}" if 'species' in rank_boot else '')
            writer.writerow(row)


if __name__ == '__main__':
    main()
