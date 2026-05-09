#!/usr/bin/env python3
"""
Build a Nextflow samplesheet CSV from a directory of FASTQ files.

Expected filename format (Illumina convention):
    {sample_id}-{marker}_{S-index}_{R1|R2}_001.fastq.gz

Examples:
    H4D-EA100_2024-18S_S244_R1_001.fastq.gz
    K9G-EA71_2024-rbcL_S503_R1_001.fastq.gz
    P8B-Blank2-SW-rbcL_S762_R2_001.fastq.gz

Outputs:
    samplesheet.csv   -- sample, fastq_1, fastq_2, marker  (pass to --input)
    metadata.tsv      -- filtered metadata aligned to found samples
                         (only written when --metadata is supplied;
                          pass to --metadata for ecology analysis)
"""

import argparse
import csv
import re
import sys
from pathlib import Path

# Canonical marker names keyed by their upper-case representation in filenames
MARKER_ALIASES = {
    '16S':  '16S',
    '18S':  '18S',
    'ITS':  'ITS',
    'ITS1': 'ITS',
    'ITS2': 'ITS',
    'CO1':  'CO1',
    'COI':  'CO1',
    '12S':  '12S',
    'RBCL': 'RBCL',
}

# Matches the Illumina suffix: _S<number>_R1/R2_001.fastq.gz
ILLUMINA_SUFFIX = re.compile(r'_(S\d+)_(R[12])_001\.fastq\.gz$', re.IGNORECASE)


def detect_marker(stem: str):
    """Return (sample_id, canonical_marker) by finding a known marker token.

    Searches last '-' then last '_' separated token in *stem* (the filename
    after the Illumina suffix has been stripped).
    """
    for sep in ('-', '_'):
        idx = stem.rfind(sep)
        if idx == -1:
            continue
        candidate = stem[idx + 1:].upper()
        if candidate in MARKER_ALIASES:
            return stem[:idx], MARKER_ALIASES[candidate]
    return None, None


def parse_fastq_dir(fastq_dir: Path, absolute: bool):
    """Scan directory; return dict keyed by (sample_id, marker) → {R1: path, R2: path}."""
    samples: dict = {}
    skipped: list = []

    for f in sorted(fastq_dir.glob('*.fastq.gz')):
        m = ILLUMINA_SUFFIX.search(f.name)
        if not m:
            skipped.append(f.name)
            continue

        read = m.group(2).upper()
        stem = f.name[: m.start()]
        sample_id, marker = detect_marker(stem)

        if sample_id is None:
            skipped.append(f.name)
            continue

        key = (sample_id, marker)
        if key not in samples:
            samples[key] = {}
        if read in samples[key]:
            sys.exit(
                f"ERROR: duplicate {read} for sample='{sample_id}' marker='{marker}'\n"
                f"  existing: {samples[key][read]}\n"
                f"  new:      {f}"
            )
        samples[key][read] = f.resolve() if absolute else f

    return samples, skipped


def load_metadata(path: Path):
    """Return (header list, dict keyed by sample_id) from a TSV/CSV metadata file."""
    delimiter = '\t' if path.suffix.lower() in ('.tsv', '.txt') else ','
    rows = {}
    header = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        header = reader.fieldnames or []
        # Accept 'sample_id', 'sample', or 'SampleID' as the key column
        id_col = next(
            (c for c in header if c.lower() in ('sample_id', 'sample', 'sampleid')),
            header[0] if header else None
        )
        if id_col is None:
            sys.exit("ERROR: cannot identify sample-ID column in metadata file")
        for row in reader:
            rows[row[id_col]] = row
    return id_col, header, rows


def write_samplesheet(samples, output_path: Path):
    with output_path.open('w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['sample', 'fastq_1', 'fastq_2', 'marker'])
        for (sample_id, marker), reads in sorted(samples.items()):
            r1 = reads.get('R1', '')
            r2 = reads.get('R2', '')
            writer.writerow([sample_id, r1, r2, marker])


def write_filtered_metadata(samples, id_col, header, meta_rows, output_path: Path):
    """Write metadata rows whose sample_id matches a parsed sample; warn on gaps."""
    sample_ids = {sid for sid, _ in samples}
    matched = {sid for sid in sample_ids if sid in meta_rows}
    no_meta = sample_ids - matched
    no_fastq = set(meta_rows) - sample_ids

    if no_meta:
        print(
            f"WARNING: {len(no_meta)} sample(s) have FASTQ files but no metadata row:\n"
            + ''.join(f'  {s}\n' for s in sorted(no_meta)),
            file=sys.stderr,
        )
    if no_fastq:
        print(
            f"WARNING: {len(no_fastq)} metadata row(s) have no matching FASTQ files:\n"
            + ''.join(f'  {s}\n' for s in sorted(no_fastq)),
            file=sys.stderr,
        )

    delimiter = '\t' if output_path.suffix.lower() in ('.tsv', '.txt') else ','
    with output_path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter=delimiter,
                                extrasaction='ignore')
        writer.writeheader()
        for sid in sorted(matched):
            writer.writerow(meta_rows[sid])

    return len(matched)


def main():
    parser = argparse.ArgumentParser(
        description='Create a Nextflow eDNA samplesheet from a FASTQ directory',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('fastq_dir', help='Directory containing *.fastq.gz files')
    parser.add_argument(
        '-o', '--output', default='samplesheet.csv',
        help='Output samplesheet CSV path (default: samplesheet.csv)',
    )
    parser.add_argument(
        '-m', '--metadata',
        help='Metadata TSV/CSV file to merge; produces a filtered copy alongside the samplesheet',
    )
    parser.add_argument(
        '--metadata-out', default=None,
        help='Path for filtered metadata output (default: <output_stem>_metadata.tsv)',
    )
    parser.add_argument(
        '--absolute', action='store_true',
        help='Write absolute FASTQ paths in the samplesheet',
    )
    args = parser.parse_args()

    fastq_dir = Path(args.fastq_dir).resolve()
    if not fastq_dir.is_dir():
        sys.exit(f"ERROR: '{fastq_dir}' is not a directory")

    samples, skipped = parse_fastq_dir(fastq_dir, args.absolute)

    if skipped:
        print(
            f"WARNING: {len(skipped)} file(s) skipped (unrecognised pattern):",
            file=sys.stderr,
        )
        for s in skipped:
            print(f'  {s}', file=sys.stderr)

    if not samples:
        sys.exit("ERROR: no samples could be parsed from the FASTQ directory")

    output_path = Path(args.output)
    write_samplesheet(samples, output_path)
    print(f"Wrote {len(samples)} sample-marker pairs to {output_path}")

    # Summarise markers found
    marker_counts: dict = {}
    for _, marker in samples:
        marker_counts[marker] = marker_counts.get(marker, 0) + 1
    for marker, count in sorted(marker_counts.items()):
        print(f"  {marker}: {count} sample(s)")

    # Optional metadata merge
    if args.metadata:
        meta_path = Path(args.metadata)
        if not meta_path.exists():
            sys.exit(f"ERROR: metadata file not found: {meta_path}")

        meta_out = Path(args.metadata_out) if args.metadata_out else \
            output_path.with_name(output_path.stem + '_metadata.tsv')

        id_col, header, meta_rows = load_metadata(meta_path)
        n = write_filtered_metadata(samples, id_col, header, meta_rows, meta_out)
        print(f"Wrote {n} matched metadata rows to {meta_out}")
        print(f"\nNextflow command:")
        print(f"  nextflow run main.nf --input {output_path} --metadata {meta_out} \\")
        print(f"    --markers '{','.join(sorted(marker_counts))}' --outdir results")
    else:
        print(
            "\nNote: no metadata file provided — ecology analyses will run in "
            "unsupervised mode (no grouping, PERMANOVA, or differential abundance).",
            file=sys.stderr,
        )
        print(f"\nNextflow command:")
        print(f"  nextflow run main.nf --input {output_path} \\")
        print(f"    --markers '{','.join(sorted(marker_counts))}' --outdir results")


if __name__ == '__main__':
    main()
