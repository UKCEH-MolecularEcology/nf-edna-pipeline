process HONEYPI_STANDARDIZE {
    tag "${marker}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/asv_taxonomy/${marker}", mode: 'copy'

    input:
    tuple val(marker), path(counts_filtered), path(consolidated_fasta), path(taxonomy_txt), path(id_map)

    output:
    tuple val(marker), path("${marker}.merged_asv_table.tsv"),     emit: merged_table
    tuple val(marker), path("${marker}.taxonomy_by_sequence.tsv"), emit: taxonomy_by_sequence
    path 'versions.yml',                                           emit: versions

    script:
    def prefix = "${marker}"
    """
    #!/usr/bin/env python3

    # ── ASV_ID -> sequence, from the consolidated FASTA ─────────────────────
    asv_to_seq = {}
    with open("${consolidated_fasta}") as fh:
        asv_id = None
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                asv_id = line[1:].split()[0]
            elif asv_id:
                asv_to_seq[asv_id] = asv_to_seq.get(asv_id, '') + line

    # ── R-sanitised column name -> original (underscore) sample ID ──────────
    id_map = {}
    with open("${id_map}") as fh:
        next(fh)  # header
        for line in fh:
            parts = line.rstrip('\\n').split('\\t')
            if len(parts) == 2:
                id_map[parts[0]] = parts[1]

    # ── Sequence-keyed merged ASV table (matches merged_asv_table.tsv's
    #    convention across every other marker: the 'asv_id' column actually
    #    holds the raw sequence, so downstream ecology joins work unchanged) ──
    with open("${counts_filtered}") as fin, open("${prefix}.merged_asv_table.tsv", 'w') as fout:
        header = fin.readline().rstrip('\\n').split('\\t')
        sample_cols = header[1:]
        orig_cols = [id_map.get(c, c) for c in sample_cols]
        fout.write('\\t'.join(['asv_id'] + orig_cols) + '\\n')

        for line in fin:
            parts = line.rstrip('\\n').split('\\t')
            asv_id = parts[0]
            seq = asv_to_seq.get(asv_id, asv_id)
            fout.write('\\t'.join([seq] + parts[1:]) + '\\n')

    # ── Taxonomy string ("k__X;p__X|Y;c__X|Y|Z;...") -> rank columns,
    #    keyed by sequence to match the merged table's join key ─────────────
    rank_letter_to_col = {
        'k': 'Kingdom', 'p': 'Phylum', 'c': 'Class',
        'o': 'Order',   'f': 'Family', 'g': 'Genus', 's': 'Species',
    }
    rank_cols = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    with open("${taxonomy_txt}") as fin, open("${prefix}.taxonomy_by_sequence.tsv", 'w') as fout:
        next(fin)  # header: ASV_ID  taxonomy  confidence
        fout.write('\\t'.join(['sequence'] + rank_cols + ['confidence']) + '\\n')

        for line in fin:
            parts = line.rstrip('\\n').split('\\t')
            if len(parts) < 2:
                continue
            asv_id, tax_string = parts[0], parts[1]
            confidence = parts[2] if len(parts) > 2 else ''
            seq = asv_to_seq.get(asv_id, asv_id)

            ranks = dict.fromkeys(rank_cols, '')
            for field in tax_string.split(';'):
                if '__' not in field:
                    continue
                letter, path_value = field.split('__', 1)
                col = rank_letter_to_col.get(letter.strip().lower())
                if col:
                    ranks[col] = path_value.split('|')[-1]

            fout.write('\\t'.join([seq] + [ranks[c] for c in rank_cols] + [confidence]) + '\\n')

    with open('versions.yml', 'w') as fh:
        fh.write('"HONEYPI_STANDARDIZE":\\n')
    """

    stub:
    def prefix = "${marker}"
    """
    touch ${prefix}.merged_asv_table.tsv
    touch ${prefix}.taxonomy_by_sequence.tsv
    printf '"HONEYPI_STANDARDIZE":\n' > versions.yml
    """
}
