process HONEYPI_MERGE_DUPLICATES {
    tag "merge_duplicates"
    label 'process_low'

    container 'python:3.10-slim'

    publishDir "${params.outdir}/pits/counts", mode: 'copy'

    input:
    path counts_txt    // ASVs_counts_filtered.txt
    path taxonomy_txt  // ASVs_taxonomy.txt

    output:
    path 'ASVs_counts_merged.txt', emit: counts
    path 'versions.yml',           emit: versions

    script:
    """
    #!/usr/bin/env python3
    from collections import defaultdict

    # Load taxonomy: ASV_ID -> taxonomy_string
    tax_map = {}
    with open("${taxonomy_txt}") as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip('\\n').split('\\t')
            if len(parts) >= 2:
                tax_map[parts[0]] = parts[1]

    # Load count table
    samples    = None
    tax_counts = defaultdict(lambda: None)

    with open("${counts_txt}") as fh:
        col_header = fh.readline().rstrip('\\n').split('\\t')
        # First column is ASV_ID, rest are sample names
        samples = col_header[1:]

        for line in fh:
            parts  = line.rstrip('\\n').split('\\t')
            asv_id = parts[0]
            counts = [int(x) if x.isdigit() else 0 for x in parts[1:]]

            taxonomy = tax_map.get(asv_id, 'k__unclassified')
            if tax_counts[taxonomy] is None:
                tax_counts[taxonomy] = [0] * len(samples)
            for i, c in enumerate(counts):
                tax_counts[taxonomy][i] += c

    # Write merged table sorted by total abundance (descending)
    rows = [(t, c) for t, c in tax_counts.items() if c is not None]
    rows.sort(key=lambda x: -sum(x[1]))

    with open('ASVs_counts_merged.txt', 'w') as fout:
        fout.write('\\t'.join(['taxonomy'] + samples) + '\\n')
        for taxonomy, counts in rows:
            fout.write('\\t'.join([taxonomy] + [str(c) for c in counts]) + '\\n')

    n_input  = len(tax_counts)
    n_unique = len(rows)
    print(f"Merged {n_input} taxa groups -> {n_unique} unique taxonomy entries", flush=True)

    with open('versions.yml', 'w') as fh:
        fh.write('"HONEYPI_MERGE_DUPLICATES":\\n')
    """

    stub:
    """
    printf 'taxonomy\tsample-A\tsample-B\nk__Viridiplantae;p__Magnoliophyta;c__Liliopsida;o__Poales;f__Poaceae;g__Hordeum\t100\t80\nk__unclassified\t50\t120\n' > ASVs_counts_merged.txt
    printf '"HONEYPI_MERGE_DUPLICATES":\n' > versions.yml
    """
}
