process TAPIRS_BLAST_LCA {
    tag "${meta.id}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    input:
    tuple val(meta), path(blast_tsv)
    path taxdump_dir      // NCBI taxdump dir (nodes/names/merged.dmp), or assets/NO_FILE when lineage_mode=='coidb'
    val lineage_mode      // 'blast' (NCBI staxids via taxdump) | 'coidb' (self-describing reference headers)
    val mlca_opts         // { identity, coverage, bitscore, majority, min_hits }

    output:
    tuple val(meta), path('*.lca.tsv'), emit: lca
    path 'versions.yml',                emit: versions

    script:
    def prefix = "${meta.id}"
    def m = mlca_opts
    def taxdump_arg = lineage_mode == 'coidb' ? '' : "--taxdump ${taxdump_dir} "
    """
    taxdump_lineage.py \\
        ${taxdump_arg}\\
        --mode ${lineage_mode} \\
        --in ${blast_tsv} \\
        --out ${prefix}.blast.tax.tsv

    mlca.py \\
        --in ${prefix}.blast.tax.tsv \\
        --out ${prefix}.lca.tsv \\
        --identity ${m.identity} \\
        --coverage ${m.coverage} \\
        --bitscore ${m.bitscore} \\
        --majority ${m.majority} \\
        --min-hits ${m.min_hits}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.lca.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
