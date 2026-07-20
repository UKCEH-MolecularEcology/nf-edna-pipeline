process TAPIRS_BLAST {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/blast:2.14.1--pl5321h6f7f691_0'

    input:
    tuple val(meta), path(nonchimeras)
    path blast_db_dir      // directory containing the BLAST db files
    val blast_db_prefix    // db basename (files under blast_db_dir named <prefix>.n**)
    val blast_opts         // { min_perc_ident, min_evalue, max_target_seqs }

    output:
    tuple val(meta), path('*.blast.tsv'), emit: blast_tsv
    path 'versions.yml',                  emit: versions

    script:
    def prefix = "${meta.id}"
    def b = blast_opts
    // blastn's own default -task is 'megablast' (tuned for near-identical
    // matches), not the more sensitive 'blastn' task -- for divergent
    // cross-species barcode matching (e.g. CO1 vs a broad reference like
    // nt) megablast finds essentially nothing. Defaults to megablast here
    // to keep 12S's already-working behavior unchanged; set blast.task in
    // a marker's config to override.
    def blast_task = b.task ?: 'megablast'
    """
    if [ -s ${nonchimeras} ]; then
        blastn -query ${nonchimeras} \\
            -db ${blast_db_dir}/${blast_db_prefix} \\
            -task ${blast_task} \\
            -outfmt "6 qseqid stitle sacc staxids pident qcovs evalue bitscore" \\
            -perc_identity ${b.min_perc_ident} \\
            -evalue ${b.min_evalue} \\
            -max_target_seqs ${b.max_target_seqs} \\
            -num_threads ${task.cpus} \\
            -out ${prefix}.blast.tsv
    else
        touch ${prefix}.blast.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -1 | sed 's/^blastn: //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.blast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.14.1
    END_VERSIONS
    """
}
