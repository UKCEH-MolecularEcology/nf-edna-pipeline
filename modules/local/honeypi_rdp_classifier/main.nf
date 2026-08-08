process HONEYPI_RDP_CLASSIFIER {
    tag "rdp_classifier"
    label 'process_medium'

    container 'quay.io/biocontainers/rdp_classifier:2.13--hdfd78af_1'

    publishDir "${params.outdir}/pits/rdp_classifier", mode: 'copy'

    input:
    path  fasta
    path  db_dir
    val   confidence_threshold

    output:
    path 'assigned_taxonomy_prelim.txt', emit: prelim
    path 'ASVs_taxonomy.txt',            emit: taxonomy
    path 'versions.yml',                 emit: versions

    script:
    def conf = confidence_threshold ?: 0.5
    """
    db_properties=\$(find -L "${db_dir}" -name '*.properties' | head -1)
    if [ -z "\${db_properties}" ]; then
        echo "ERROR: No .properties file found in ${db_dir}" >&2; exit 1
    fi

    n_seqs=\$(grep -c '^>' "${fasta}" || echo 0)
    echo "Classifying \${n_seqs} sequences against \${db_properties}"

    rdp_classifier classify \\
        -Xmx${task.memory.toGiga()}g \\
        -t "\${db_properties}" \\
        -o assigned_taxonomy_prelim.txt \\
        "${fasta}"

    # Reformat RDP output to TSV taxonomy table (awk only — no python3 in this container)
    awk -v threshold=${conf} -f "${projectDir}/bin/parse_rdp.awk" assigned_taxonomy_prelim.txt > ASVs_taxonomy.txt

    printf '"HONEYPI_RDP_CLASSIFIER":\\n    rdp_classifier: 2.13\\n' > versions.yml
    """

    stub:
    """
    printf 'ASV_0000000001\t+\troot\trootrank\t1.0\tViridiplantae\tdomain\t0.99\tMagnoliophyta\tphylum\t0.95\tLiliopsida\tclass\t0.80\tPoales\torder\t0.75\tPoaceae\tfamily\t0.70\tHordeum\tgenus\t0.65\n' > assigned_taxonomy_prelim.txt
    printf 'ASV_ID\ttaxonomy\tconfidence\nASV_0000000001\tk__Viridiplantae;p__Magnoliophyta;c__Liliopsida;o__Poales;f__Poaceae;g__Hordeum\t0.65\nASV_0000000002\tk__unclassified\t1.0\n' > ASVs_taxonomy.txt
    printf '"HONEYPI_RDP_CLASSIFIER":\n    rdp_classifier: 2.13\n' > versions.yml
    """
}
