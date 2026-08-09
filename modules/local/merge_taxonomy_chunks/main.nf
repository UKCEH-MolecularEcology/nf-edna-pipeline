process MERGE_TAXONOMY_CHUNKS {
    tag "${marker}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/taxonomy/${marker}", mode: 'copy'

    input:
    val  marker
    path chunks   // one *.taxonomy.tsv per FASTA chunk from TAXONOMY (dada2 method)

    output:
    tuple val(marker), path("${marker}.taxonomy.tsv"), emit: taxonomy

    script:
    """
    awk 'FNR==1 && NR!=1 { next } { print }' ${chunks} > ${marker}.taxonomy.tsv
    """

    stub:
    """
    touch ${marker}.taxonomy.tsv
    """
}
