process TAXONOMY {
    tag "${meta.id}_${marker}"
    label 'process_high'

    container {
        tax_method == 'blast'
            ? 'biocontainers/blast:2.14.1--pl5321h6f7f691_0'
            : tax_method == 'rdp'
                ? 'quay.io/biocontainers/rdp-classifier:2.13--hdfd78af_1'
                : 'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'
    }

    publishDir "${params.outdir}/taxonomy/${marker}", mode: 'copy'

    input:
    tuple val(meta), path(fasta)
    path  tax_db
    val   db_type    // silva | unite | midori | pr2 | rdp
    val   tax_method // dada2 | blast | rdp
    val   marker

    output:
    tuple val(meta), path('*.taxonomy.tsv'),  emit: taxonomy
    path 'versions.yml',                       emit: versions

    script:
    def prefix = "${meta.id}_${marker}"

    if (tax_method == 'rdp') {
        def rdp_xmx  = params.rdp_mem       ?: '8g'
        def rdp_boot = params.rdp_bootstrap  ?: 0.70

        """
        rdp_classifier -Xmx${rdp_xmx} classify \\
            -t ${tax_db}/rRNAClassifier.properties \\
            -f fixrank \\
            -o ${prefix}_rdp_raw.txt \\
            ${fasta}

        parse_rdp.py ${prefix}_rdp_raw.txt ${prefix}.taxonomy.tsv ${rdp_boot}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rdp_classifier: \$(rdp_classifier --version 2>&1 | grep -oP '\\d+\\.\\d+' | head -1 || echo "2.13")
            python: \$(python3 --version | cut -d' ' -f2)
        END_VERSIONS
        """

    } else if (tax_method == 'blast') {
        """
        # BLAST against curated database (pre-formatted)
        blastn \\
            -query ${fasta} \\
            -db ${tax_db} \\
            -out ${prefix}.blast_raw.tsv \\
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames" \\
            -max_target_seqs 5 \\
            -perc_identity 80 \\
            -num_threads ${task.cpus}

        # Take top hit per ASV and format as taxonomy table
        awk 'BEGIN{OFS="\\t"; print "asv_id","subject","pident","evalue","bitscore","taxonomy"}
             !seen[\$1]++{print \$1,\$2,\$3,\$11,\$12,\$14}' ${prefix}.blast_raw.tsv \\
        > ${prefix}.taxonomy.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blast: \$(blastn -version | head -1 | sed 's/blastn: //')
        END_VERSIONS
        """
    } else {
        // DADA2 assignTaxonomy (naive Bayesian classifier)
        // Handles SILVA, UNITE, MIDORI, PR2 formatted databases
        def species_db = db_type == 'silva'
            ? file(tax_db.toString().replaceAll('train_set', 'species_assignment'))
            : null
        def add_species = (db_type == 'silva' && species_db?.exists()) ? 'TRUE' : 'FALSE'

        """
        #!/usr/bin/env Rscript
        library(dada2)

        marker     <- "${marker}"
        db_type    <- "${db_type}"
        fasta_file <- "${fasta}"
        tax_db     <- "${tax_db}"
        prefix     <- "${prefix}"
        add_species <- as.logical("${add_species}")

        seqs <- getSequences(fasta_file)
        names(seqs) <- gsub(">", "", system(paste("grep '^>' ", fasta_file), intern=TRUE))
        names(seqs) <- sub(" .*", "", names(seqs))

        # Assign taxonomy
        taxa <- assignTaxonomy(
            seqs,
            tax_db,
            multithread  = ${task.cpus},
            minBoot      = 50,
            outputBootstraps = TRUE,
            tryRC        = TRUE
        )

        tax_table   <- taxa\$tax
        boot_table  <- taxa\$boot

        # Format ranks based on database type
        if (db_type %in% c("silva", "midori")) {
            rank_names <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
        } else if (db_type == "pr2") {
            rank_names <- c("Domain","Supergroup","Division","Subdivision",
                           "Class","Order","Family","Genus","Species")
        } else if (db_type == "unite") {
            rank_names <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
        } else {
            rank_names <- colnames(tax_table)
        }
        colnames(tax_table) <- rank_names[seq_len(ncol(tax_table))]

        out <- data.frame(
            asv_id = rownames(tax_table),
            tax_table,
            stringsAsFactors = FALSE
        )

        # Add confidence scores for Genus and Species
        if ("Genus" %in% colnames(boot_table)) {
            out\$genus_boot   <- boot_table[,"Genus"]
        }
        if ("Species" %in% colnames(boot_table)) {
            out\$species_boot <- boot_table[,"Species"]
        }

        write.table(out, paste0(prefix, ".taxonomy.tsv"),
                    sep = "\\t", quote = FALSE, row.names = FALSE, na = "")

        writeLines(
            c(
                paste0('"${task.process}":'),
                paste0('    dada2: ', packageVersion('dada2')),
                paste0('    R: ', R.version\$major, '.', R.version\$minor)
            ),
            "versions.yml"
        )
        """
    }

    stub:
    def prefix = "${meta.id}_${meta.marker}"
    """
    touch ${prefix}.taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.3
    END_VERSIONS
    """

}
