process TAPIRS_OTU_TABLE {
    tag "${out_prefix}"
    label 'process_medium'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/tapirs", mode: 'copy'

    input:
    path tax_files      // collected across ALL samples for this marker: .lca.tsv or .krk.tax.tsv
    path rerep_files     // collected across ALL samples for this marker: .rerep.fasta
    val out_prefix
    val lowest_rank
    val highest_rank

    output:
    path "${out_prefix}.tsv",              emit: otu_table
    path "${out_prefix}_full_lineage.tsv", emit: otu_table_full_lineage
    path 'versions.yml',                   emit: versions

    script:
    """
    tapirs_otu_table.py \\
        --tax-files ${tax_files} \\
        --rerep-files ${rerep_files} \\
        --lowest-rank ${lowest_rank} \\
        --highest-rank ${highest_rank} \\
        --out-regular ${out_prefix}.tsv \\
        --out-full-lineage ${out_prefix}_full_lineage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${out_prefix}.tsv
    touch ${out_prefix}_full_lineage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
