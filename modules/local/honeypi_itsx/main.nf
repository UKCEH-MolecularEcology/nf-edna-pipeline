process HONEYPI_ITSX {
    tag "itsx_${its_region}"
    label 'process_medium'

    container 'quay.io/biocontainers/itsx:1.1.3--hdfd78af_1'

    publishDir "${params.outdir}/pits/itsx", mode: 'copy'

    input:
    path asvs_fasta
    val  its_region
    val  its_min_len
    val  its_max_len

    output:
    path 'ASVs_ITSxed.fasta',  emit: itsxed
    path 'itsx_out.summary.txt', emit: summary
    path 'versions.yml',       emit: versions

    script:
    """
    # Count input sequences
    n_seqs=\$(grep -c '^>' "${asvs_fasta}" || echo 0)
    echo "ITSx input: \${n_seqs} sequences"

    ITSx \\
        -i "${asvs_fasta}" \\
        -o itsx_out \\
        --taxa all \\
        --complement T \\
        -t ${its_region} \\
        --cpu ${task.cpus} \\
        --graphical F \\
        --save_regions ${its_region} \\
        --partial 50 || true

    # Collect the extracted region
    ITS_FILE="itsx_out.${its_region}.fasta"
    if [ ! -f "\${ITS_FILE}" ] || [ ! -s "\${ITS_FILE}" ]; then
        echo "WARNING: ITSx produced no ${its_region} sequences. Creating empty file."
        touch ASVs_ITSxed.fasta
    else
        # Length filter using Python (replaces vsearch)
        python3 << 'PYEOF'
import sys

min_len = ${its_min_len}
max_len = ${its_max_len}
its_file = "itsx_out.${its_region}.fasta"

kept = 0
removed = 0
current_hdr = None
current_seq  = []

def write_record(out, hdr, seq):
    global kept, removed
    s = ''.join(seq)
    if min_len <= len(s) <= max_len:
        out.write(hdr + '\\n' + s + '\\n')
        kept += 1
    else:
        removed += 1

with open(its_file) as fin, open('ASVs_ITSxed.fasta', 'w') as fout:
    for line in fin:
        line = line.rstrip()
        if line.startswith('>'):
            if current_hdr is not None:
                write_record(fout, current_hdr, current_seq)
            current_hdr = line
            current_seq  = []
        else:
            current_seq.append(line)
    if current_hdr is not None:
        write_record(fout, current_hdr, current_seq)

print(f"Length filter ({min_len}-{max_len} bp): kept {kept}, removed {removed}", flush=True)
PYEOF
    fi

    n_out=\$(grep -c '^>' ASVs_ITSxed.fasta 2>/dev/null || echo 0)
    echo "ITSx output: \${n_out} sequences after length filter"

    [ -f itsx_out.summary.txt ] || touch itsx_out.summary.txt

    ITSx --version 2>&1 | head -1 | sed 's/ITSx//' | tr -d 'v' | \\
        awk '{print "\\"HONEYPI_ITSX\\":\\n    ITSx: " \$1}' > versions.yml || \\
        printf '"HONEYPI_ITSX":\\n    ITSx: 1.1.3\\n' > versions.yml
    """

    stub:
    """
    printf '>ASV_0000000001\nACGTACGT\n' > ASVs_ITSxed.fasta
    touch itsx_out.summary.txt
    printf '"HONEYPI_ITSX":\n    ITSx: 1.1.3\n' > versions.yml
    """
}
