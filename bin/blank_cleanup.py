#!/usr/bin/env python3
"""
LOD (limit of detection) blank cleanup for the collapsed 12S SINTAX/CLARE
abundance matrix. Ported from 12S-edna-dada2-tapirs-workflow's
06_blank_cleanup_from_workflow.R.

Classifies each SAMPLE COLUMN NAME (not a metadata file) into
pcr_blank / extraction_blank / site_blank / sample via priority-ordered
regexes, computes lab-blank and site-blank LOD (mean + N*SD) and LOQ
(mean + M*SD) per taxon (+ site for site blanks), and zeroes any real-sample
cell below its taxon's LOD in three variants (lab-only, site-only, both).

Usage:
    blank_cleanup.py --clare-abundance asv_taxonomy_abundance_CLARE.csv \\
        --ncl-matrix ncl_matrix_raw.csv \\
        --pcr-blank-regex 'PCR_BLANK|PCRBLANK|NEG|POS' \\
        --extraction-blank-regex 'EXTRACTION_BLANK|EXTRACTIONBLANK|_EB_|_EB[0-9]' \\
        --site-blank-regex 'BLANK' \\
        --lod-sd-multiplier 3 --loq-sd-multiplier 10 \\
        [--enable-taxon-exclusion] [--excluded-taxa Homo_sapiens,Sus_scrofa,...]
"""

import argparse
import re
import sys

import pandas as pd

RANK_PRIORITY = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]


def classify_sample_type(name, pcr_re, extraction_re, site_re):
    nm = name.upper()
    if re.search(pcr_re, nm):
        return "pcr_blank"
    if re.search(extraction_re, nm):
        return "extraction_blank"
    if re.search(site_re, nm):
        return "site_blank"
    return "sample"


def get_site_code(name):
    nm = name.upper()
    m = re.search(r"_([A-Z]{4})BLANK", nm)
    if m:
        return m.group(1)
    m = re.search(r"_R[12]_([A-Z]{4})\d", nm)
    if m:
        return m.group(1)
    return None


def lowest_rank_conf(row):
    for rank in RANK_PRIORITY:
        val = row.get(rank)
        if val and str(val).strip():
            return rank, val, row.get(f"{rank}_conf")
    return None, None, None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clare-abundance", required=True)
    parser.add_argument("--ncl-matrix", required=True)
    parser.add_argument("--pcr-blank-regex", required=True)
    parser.add_argument("--extraction-blank-regex", required=True)
    parser.add_argument("--site-blank-regex", required=True)
    parser.add_argument("--lod-sd-multiplier", type=float, default=3)
    parser.add_argument("--loq-sd-multiplier", type=float, default=10)
    parser.add_argument("--enable-taxon-exclusion", action="store_true")
    parser.add_argument("--excluded-taxa", default="")
    args = parser.parse_args()

    clare = pd.read_csv(args.clare_abundance, dtype=str)
    ncl = pd.read_csv(args.ncl_matrix)

    sample_cols = [c for c in ncl.columns if c not in ("species", "taxonomy")]

    # Step 2: classify each sample column
    sample_info = pd.DataFrame({"sample": sample_cols})
    sample_info["sample_type"] = sample_info["sample"].apply(
        lambda s: classify_sample_type(s, args.pcr_blank_regex, args.extraction_blank_regex, args.site_blank_regex)
    )
    sample_info["site"] = sample_info["sample"].apply(get_site_code)
    sample_info.to_csv("cleanup_sample_classification.csv", index=False)

    # Step 3: lowest resolved rank/name/conf per unique taxonomy string (from CLARE)
    lookup_rows = {}
    for _, row in clare.iterrows():
        taxonomy = row.get("taxonomy", "")
        if taxonomy in lookup_rows:
            continue
        rank, name, conf = lowest_rank_conf(row)
        lookup_rows[taxonomy] = {"lowest_rank": rank, "lowest_name": name, "lowest_conf": conf}
    tax_lookup = pd.DataFrame(
        [{"taxonomy": k, **v} for k, v in lookup_rows.items()]
    )
    ncl_tax = ncl[["species", "taxonomy"]].merge(tax_lookup, on="taxonomy", how="left")
    ncl_tax.to_csv("cleanup_taxonomy_lookup.csv", index=False)

    # Step 4: long format
    long_df = ncl.melt(id_vars=["species", "taxonomy"], value_vars=sample_cols,
                        var_name="sample", value_name="reads")
    long_df = long_df.merge(sample_info, on="sample", how="left")
    long_df = long_df.merge(ncl_tax[["species", "taxonomy", "lowest_rank", "lowest_name", "lowest_conf"]],
                             on=["species", "taxonomy"], how="left")

    # Step 5: blank stats / LOD-LOQ
    lab_blanks = long_df[long_df["sample_type"].isin(["pcr_blank", "extraction_blank"])]
    lab_stats = lab_blanks.groupby("species")["reads"].agg(
        mean_lab_blank_reads="mean", sd_lab_blank_reads="std", max_lab_blank_reads="max"
    ).reset_index()
    lab_stats["sd_lab_blank_reads"] = lab_stats["sd_lab_blank_reads"].fillna(0)
    lab_stats["lab_LOD"] = lab_stats["mean_lab_blank_reads"] + args.lod_sd_multiplier * lab_stats["sd_lab_blank_reads"]
    lab_stats["lab_LOQ"] = lab_stats["mean_lab_blank_reads"] + args.loq_sd_multiplier * lab_stats["sd_lab_blank_reads"]
    lab_stats.to_csv("cleanup_lab_blank_stats.csv", index=False)

    site_blanks = long_df[(long_df["sample_type"] == "site_blank") & long_df["site"].notna()]
    site_stats = site_blanks.groupby(["species", "site"])["reads"].agg(
        mean_site_blank_reads="mean", sd_site_blank_reads="std", max_site_blank_reads="max"
    ).reset_index()
    site_stats["sd_site_blank_reads"] = site_stats["sd_site_blank_reads"].fillna(0)
    site_stats["site_LOD"] = site_stats["mean_site_blank_reads"] + args.lod_sd_multiplier * site_stats["sd_site_blank_reads"]
    site_stats["site_LOQ"] = site_stats["mean_site_blank_reads"] + args.loq_sd_multiplier * site_stats["sd_site_blank_reads"]
    site_stats.to_csv("cleanup_site_blank_stats.csv", index=False)

    # Step 6: diagnostic table (real samples only)
    real = long_df[long_df["sample_type"] == "sample"]
    diag = real.groupby(["species", "taxonomy", "lowest_rank", "lowest_name", "lowest_conf", "site"], dropna=False)["reads"].agg(
        max_sample_reads="max", mean_sample_reads="mean", n_positive_samples=lambda s: int((s > 0).sum())
    ).reset_index()
    diag = diag.merge(lab_stats[["species", "max_lab_blank_reads"]], on="species", how="left")
    diag = diag.merge(site_stats[["species", "site", "max_site_blank_reads"]], on=["species", "site"], how="left")
    diag["lab_blank_to_sample_ratio"] = diag.apply(
        lambda r: (r["max_lab_blank_reads"] / r["max_sample_reads"]) if r["max_sample_reads"] else pd.NA, axis=1)
    diag["site_blank_to_sample_ratio"] = diag.apply(
        lambda r: (r["max_site_blank_reads"] / r["max_sample_reads"]) if r["max_sample_reads"] else pd.NA, axis=1)
    diag = diag.sort_values(
        by=["lab_blank_to_sample_ratio", "site_blank_to_sample_ratio"], ascending=False, na_position="last"
    )
    diag.to_csv("cleanup_blank_diagnostic_table.csv", index=False)

    # Step 7: apply LOD cleanup
    sample_long = real.merge(lab_stats[["species", "lab_LOD"]], on="species", how="left")
    sample_long = sample_long.merge(site_stats[["species", "site", "site_LOD"]], on=["species", "site"], how="left")
    sample_long["lab_LOD"] = sample_long["lab_LOD"].fillna(0)
    sample_long["site_LOD"] = sample_long["site_LOD"].fillna(0)
    sample_long["reads_lab_clean"] = sample_long.apply(
        lambda r: 0 if r["reads"] < r["lab_LOD"] else r["reads"], axis=1)
    sample_long["reads_site_clean"] = sample_long.apply(
        lambda r: 0 if r["reads"] < r["site_LOD"] else r["reads"], axis=1)
    sample_long["reads_both_clean"] = sample_long.apply(
        lambda r: 0 if r["reads"] < max(r["lab_LOD"], r["site_LOD"]) else r["reads"], axis=1)

    # Step 8: human/livestock exclusion -- OFF by default, matching the
    # current Snakemake runtime behavior (present in that pipeline's code
    # but commented out there too).
    if args.enable_taxon_exclusion and args.excluded_taxa:
        excluded = args.excluded_taxa.split(",")
        pattern = "|".join(re.escape(t) for t in excluded)
        sample_long = sample_long[~sample_long["taxonomy"].fillna("").str.contains(pattern, case=False, regex=True)]

    # Step 9: wide cleaned matrices
    def make_matrix(value_col, out_path):
        wide = sample_long.pivot_table(
            index=["species", "taxonomy", "lowest_rank", "lowest_name", "lowest_conf"],
            columns="sample", values=value_col, fill_value=0, aggfunc="sum"
        ).reset_index()
        wide.to_csv(out_path, sep="\t", index=False)

    make_matrix("reads_lab_clean", "ncl_cleaned_labLOD.csv")
    make_matrix("reads_site_clean", "ncl_cleaned_siteLOD.csv")
    make_matrix("reads_both_clean", "ncl_cleaned_bothLOD.csv")

    # Step 10: presence/absence versions
    def make_pa(value_col, out_path):
        pa = sample_long.copy()
        pa[value_col] = (pa[value_col] > 0).astype(int)
        wide = pa.pivot_table(
            index=["species", "taxonomy", "lowest_rank", "lowest_name", "lowest_conf"],
            columns="sample", values=value_col, fill_value=0, aggfunc="max"
        ).reset_index()
        wide.to_csv(out_path, sep="\t", index=False)

    make_pa("reads_lab_clean", "ncl_cleaned_labLOD_pa.csv")
    make_pa("reads_site_clean", "ncl_cleaned_siteLOD_pa.csv")
    make_pa("reads_both_clean", "ncl_cleaned_bothLOD_pa.csv")

    # Step 11: long output
    long_out_cols = ["species", "taxonomy", "lowest_rank", "lowest_name", "lowest_conf",
                      "sample", "site", "reads", "reads_lab_clean", "reads_site_clean", "reads_both_clean"]
    sample_long[long_out_cols].to_csv("ncl_cleaned_long.csv", sep="\t", index=False)

    # Step 12: summary
    summary = pd.DataFrame([{
        "total_raw_reads": sample_long["reads"].sum(),
        "total_lab_clean_reads": sample_long["reads_lab_clean"].sum(),
        "total_site_clean_reads": sample_long["reads_site_clean"].sum(),
        "total_both_clean_reads": sample_long["reads_both_clean"].sum(),
        "n_nonzero_raw": int((sample_long["reads"] > 0).sum()),
        "n_nonzero_lab_clean": int((sample_long["reads_lab_clean"] > 0).sum()),
        "n_nonzero_site_clean": int((sample_long["reads_site_clean"] > 0).sum()),
        "n_nonzero_both_clean": int((sample_long["reads_both_clean"] > 0).sum()),
    }])
    summary.to_csv("cleanup_summary.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())
