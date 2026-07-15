#!/usr/bin/env python3
"""
Merge per-database SINTAX parsed-taxonomy tables into one wide comparison
table plus a per-database summary. Ported from
12S-edna-dada2-tapirs-workflow's 03b_sintax_merge.R.

Usage:
    sintax_merge.py --asv-lookup asv_lookup.tsv \\
        --parsed INBO=inbo_parsed.tsv MIDORI=midori_parsed.tsv CLARE=clare_parsed.tsv \\
        --out-compare asv_taxonomy_compare.tsv --out-summary sintax_database_summary.tsv
"""

import argparse
import sys

import pandas as pd


def summarise(df, db_name):
    n_asv = len(df)

    def frac_or_zero(series):
        return series.mean(skipna=True) if len(series.dropna()) else float("nan")

    return {
        "database": db_name,
        "n_asv": n_asv,
        "assigned_family": int((df["family"].fillna("") != "").sum()),
        "assigned_genus": int((df["genus"].fillna("") != "").sum()),
        "assigned_species": int((df["species"].fillna("") != "").sum()),
        "mean_genus_conf": frac_or_zero(pd.to_numeric(df["genus_conf"], errors="coerce")),
        "mean_species_conf": frac_or_zero(pd.to_numeric(df["species_conf"], errors="coerce")),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asv-lookup", required=True)
    parser.add_argument("--parsed", nargs="+", required=True, help="DB=path pairs, e.g. INBO=inbo_parsed.tsv")
    parser.add_argument("--out-compare", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    db_paths = dict(p.split("=", 1) for p in args.parsed)

    lookup = pd.read_csv(args.asv_lookup, sep="\t", dtype=str)
    compare = lookup.copy()
    summary_rows = []

    for db_name, path in db_paths.items():
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        summary_rows.append(summarise(df, db_name))

        renamed = df[["ASV", "taxonomy", "family", "genus", "species", "genus_conf", "species_conf"]].rename(
            columns={
                "taxonomy": f"taxonomy_{db_name}",
                "family": f"family_{db_name}",
                "genus": f"genus_{db_name}",
                "species": f"species_{db_name}",
                "genus_conf": f"genus_conf_{db_name}",
                "species_conf": f"species_conf_{db_name}",
            }
        )
        compare = compare.merge(renamed, on="ASV", how="left")

    compare.to_csv(args.out_compare, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(args.out_summary, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())
