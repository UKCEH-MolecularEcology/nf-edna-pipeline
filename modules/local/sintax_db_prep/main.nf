process SINTAX_DB_PREP {
    tag "${db_name}"
    label 'process_single'

    container 'quay.io/biocontainers/vsearch:2.27.0--h6a68c12_0'

    input:
    tuple val(db_name), path(raw_fasta)

    output:
    tuple val(db_name), path('*.fasta'), emit: fasta

    script:
    def out = "${db_name}.sintax_ref.fasta"
    if (raw_fasta.name.endsWith('.gz')) {
        """
        gunzip -c ${raw_fasta} > ${out}
        """
    } else {
        """
        ln -s ${raw_fasta} ${out}
        """
    }

    stub:
    """
    touch ${db_name}.sintax_ref.fasta
    """
}
