process SINTAX_RUN {
    tag "${marker}_${db_name}"
    label 'process_medium'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/usearch12:12.0-beta1'

    input:
    tuple val(marker), path(asv_fasta), path(asv_lookup), path(merged_table)
    tuple val(db_name), path(ref_fasta)
    val cutoff

    output:
    tuple val(marker), val(db_name), path('*_parsed.tsv'),               emit: parsed
    tuple val(marker), val(db_name), path('*_asv_taxonomy.tsv'),         emit: taxonomy
    tuple val(marker), val(db_name), path('species_abundance_*.csv'),   emit: species_abundance
    tuple val(marker), val(db_name), path('asv_taxonomy_abundance_*.csv'), emit: taxonomy_abundance
    path 'versions.yml',                                                 emit: versions

    script:
    def prefix = "${marker}_${db_name}"
    """
    usearch12 -sintax ${asv_fasta} \\
        -db ${ref_fasta} \\
        -tabbedout ${prefix}.sintax.tsv \\
        -strand both \\
        -sintax_cutoff ${cutoff} \\
        -threads ${task.cpus}

    sintax_parse.py \\
        --sintax-out ${prefix}.sintax.tsv \\
        --asv-lookup ${asv_lookup} \\
        --abundance-table ${merged_table} \\
        --db-name ${db_name} \\
        --out-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch12: 12.0-beta1
    END_VERSIONS
    """

    stub:
    def prefix = "${marker}_${db_name}"
    """
    touch ${prefix}_parsed.tsv
    touch ${prefix}_asv_taxonomy.tsv
    touch species_abundance_${db_name}.csv
    touch asv_taxonomy_abundance_${db_name}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch12: 12.0-beta1
    END_VERSIONS
    """
}
