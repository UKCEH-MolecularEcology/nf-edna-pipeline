process HONEYPI_FILTER_ASV_TABLE {
    tag "filter_asv_table"
    label 'process_low'

    container 'python:3.10-slim'

    publishDir "${params.outdir}/pits/counts", mode: 'copy'

    input:
    path counts_txt   // ASVs_dada2_counts.txt (all DADA2 ASVs)
    path asvs_fasta   // consolidated ASVs.fasta (ASVs to keep)

    output:
    path 'ASVs_counts_filtered.txt', emit: counts
    path 'versions.yml',             emit: versions

    script:
    """
    #!/usr/bin/env python3

    # Collect ASV IDs from the consolidated FASTA
    keep_ids = set()
    with open("${asvs_fasta}") as fh:
        for line in fh:
            if line.startswith('>'):
                keep_ids.add(line[1:].strip().split()[0])

    print(f"Keeping {len(keep_ids)} ASV IDs from consolidated FASTA", flush=True)

    # Filter the count table
    kept = 0
    dropped = 0
    with open("${counts_txt}") as fin, open('ASVs_counts_filtered.txt', 'w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            asv_id = line.split('\\t', 1)[0].strip()
            if asv_id in keep_ids:
                fout.write(line)
                kept += 1
            else:
                dropped += 1

    print(f"Filter result: {kept} kept, {dropped} dropped", flush=True)

    with open('versions.yml', 'w') as fh:
        fh.write('"HONEYPI_FILTER_ASV_TABLE":\\n')
    """

    stub:
    """
    printf 'ASV_ID\tsample-A\tsample-B\nASV_0000000001\t100\t80\nASV_0000000002\t50\t120\n' > ASVs_counts_filtered.txt
    printf '"HONEYPI_FILTER_ASV_TABLE":\n' > versions.yml
    """
}
