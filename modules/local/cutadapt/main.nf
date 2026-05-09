process CUTADAPT {
    tag "${meta.id}_${meta.marker}"
    label 'process_medium'

    container 'biocontainers/cutadapt:4.6--py39hf95cd2a_1'

    publishDir "${params.outdir}/trimmed/${meta.marker}", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trimmed.fastq.gz'), emit: reads
    tuple val(meta), path('*.log'),              emit: log
    path 'versions.yml',                         emit: versions

    script:
    def prefix     = "${meta.id}_${meta.marker}"
    def fwd        = meta.fwd_primer ?: params.primers[meta.marker].fwd
    def rev        = meta.rev_primer ?: params.primers[meta.marker].rev
    def min_len    = meta.min_length ?: 50
    def max_len    = meta.max_length ? "--maximum-length ${meta.max_length}" : ''
    def extra_args = params.cutadapt_args ?: '--discard-untrimmed'

    if (meta.single_end) {
        """
        cutadapt \\
            -g ${fwd} \\
            --minimum-length ${min_len} \\
            ${max_len} \\
            --cores ${task.cpus} \\
            ${extra_args} \\
            -o ${prefix}.trimmed.fastq.gz \\
            ${reads[0]} \\
            > ${prefix}.log 2>&1

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        // Linked adapter strategy: trim fwd from R1, rev-complement of rev from R2
        def rev_rc = reverseComplement(rev)
        """
        cutadapt \\
            -g ${fwd} \\
            -G ${rev} \\
            -a ${rev_rc} \\
            -A ${reverseComplement(fwd)} \\
            --minimum-length ${min_len} \\
            ${max_len} \\
            --cores ${task.cpus} \\
            ${extra_args} \\
            --pair-filter=any \\
            -o ${prefix}_R1.trimmed.fastq.gz \\
            -p ${prefix}_R2.trimmed.fastq.gz \\
            ${reads[0]} ${reads[1]} \\
            > ${prefix}.log 2>&1

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }

    stub:
    def prefix = "${meta.id}_${meta.marker}"
    if (meta.single_end) {
        """
        touch ${prefix}.trimmed.fastq.gz
        touch ${prefix}.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: 4.6
        END_VERSIONS
        """
    } else {
        """
        touch ${prefix}_R1.trimmed.fastq.gz
        touch ${prefix}_R2.trimmed.fastq.gz
        touch ${prefix}.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: 4.6
        END_VERSIONS
        """
    }
}

// Groovy helper to reverse-complement a DNA primer (IUPAC)
def reverseComplement(String seq) {
    def comp = ['A':'T','T':'A','G':'C','C':'G',
                'R':'Y','Y':'R','S':'S','W':'W','K':'M','M':'K',
                'B':'V','V':'B','D':'H','H':'D','N':'N']
    seq.toUpperCase().reverse().collect { comp[it] ?: 'N' }.join('')
}
