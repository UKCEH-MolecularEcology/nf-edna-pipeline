process HONEYPI_CONSOLIDATE {
    tag "consolidate"
    label 'process_low'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/pits/consolidate", mode: 'copy'

    input:
    path dada2_asvs    // ASVs.fasta from DADA2 (all ASVs)
    path itsxed_asvs   // ASVs_ITSxed.fasta from ITSX (extracted sequences)

    output:
    path 'ASVs.fasta',   emit: fasta
    path 'versions.yml', emit: versions

    script:
    """
    #!/usr/bin/env python3

    def read_fasta(path):
        records = {}
        order   = []
        with open(path) as fh:
            hdr = None
            buf = []
            for line in fh:
                line = line.rstrip()
                if line.startswith('>'):
                    if hdr is not None:
                        records[hdr] = ''.join(buf)
                    hdr = line[1:].split()[0]
                    buf = []
                    if hdr not in records:
                        order.append(hdr)
                elif hdr:
                    buf.append(line)
            if hdr is not None:
                records[hdr] = ''.join(buf)
        return records, order

    dada2_records,  dada2_order  = read_fasta("${dada2_asvs}")
    itsxed_records, itsxed_order = read_fasta("${itsxed_asvs}")

    # Sequences NOT extracted by ITSx
    not_itsxed = {k: v for k, v in dada2_records.items() if k not in itsxed_records}

    # Combine: ITSxed first (ITS-trimmed sequences preferred), then non-ITSxed
    combined = {}
    combined.update(itsxed_records)
    combined.update(not_itsxed)

    n_itsxed     = len(itsxed_records)
    n_not_itsxed = len(not_itsxed)
    n_total      = len(combined)

    with open('ASVs.fasta', 'w') as fh:
        for asv_id in sorted(combined):
            fh.write(f'>{asv_id}\\n{combined[asv_id]}\\n')

    print(f"Consolidation: {n_itsxed} ITSxed + {n_not_itsxed} non-ITSxed = {n_total} total ASVs",
          flush=True)

    with open('versions.yml', 'w') as fh:
        fh.write('"HONEYPI_CONSOLIDATE":\\n')
    """

    stub:
    """
    printf '>ASV_0000000001\nACGTACGTACGTACGT\n>ASV_0000000002\nTGCATGCATGCATGCA\n' > ASVs.fasta
    printf '"HONEYPI_CONSOLIDATE":\n' > versions.yml
    """
}
