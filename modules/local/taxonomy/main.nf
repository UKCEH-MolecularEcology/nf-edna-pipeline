process TAXONOMY {
    tag "${marker}"
    label 'process_high'

    container {
        tax_method == 'blast'
            ? 'quay.io/biocontainers/blast:2.14.1--pl5321h6f7f691_0'
            : tax_method == 'rdp'
                ? 'quay.io/biocontainers/rdp_classifier:2.13--hdfd78af_1'
                : 'quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0'
    }

    publishDir "${params.outdir}/taxonomy/${marker}", mode: 'copy'

    input:
    tuple val(marker), path(fasta)   // marker-level, collective FASTA (post MERGE_ASV_TABLES) — see collective_taxonomy.nf
    path  tax_db
    val   db_type       // silva | unite | midori | coidb | pr2 | rdp
    val   tax_method    // dada2 | blast | rdp
    path  addspecies_db // dada2 only: exact-match species DB for addSpecies(), or assets/NO_FILE

    output:
    tuple val(marker), path('*.taxonomy.tsv'), emit: taxonomy
    path 'versions.yml',                       emit: versions

    script:
    // dada2-method calls are split into many chunks by collective_taxonomy.nf
    // (see params.tax_chunk_size) -- suffix with task.index so each chunk's
    // output filename is unique once they all land in MERGE_TAXONOMY_CHUNKS'
    // work dir together. rdp/blast run as one task on the whole merged
    // FASTA, so keep their plain marker-only filename.
    def prefix = (tax_method == 'dada2') ? "${marker}_${task.index}" : "${marker}"

    if (tax_method == 'rdp') {
        def rdp_xmx  = params.rdp_mem       ?: '8g'
        def rdp_boot = params.rdp_bootstrap  ?: 0.70

        """
        rdp_classifier -Xmx${rdp_xmx} classify \\
            -t ${tax_db}/rRNAClassifier.properties \\
            -f fixrank \\
            -o ${prefix}_rdp_raw.txt \\
            ${fasta}

        awk -v CUT="${rdp_boot}" '
        BEGIN {
            OFS="\\t"
            print "asv_id\\tKingdom\\tPhylum\\tClass\\tOrder\\tFamily\\tGenus\\tSpecies\\tgenus_boot\\tspecies_boot"
        }
        NF > 2 && length(\$1) > 0 {
            seq = \$1; split("", val); split("", boot)
            for (i = 3; i+2 <= NF; i += 3) {
                rank = tolower(\$i); val[rank] = \$(i+1); boot[rank] = \$(i+2)+0
            }
            kingdom = ""
            if (("superkingdom" in val) && boot["superkingdom"] >= CUT) kingdom = val["superkingdom"]
            else if (("kingdom" in val) && boot["kingdom"] >= CUT) kingdom = val["kingdom"]
            phylum  = (("phylum"  in val) && boot["phylum"]  >= CUT) ? val["phylum"]  : ""
            cls     = (("class"   in val) && boot["class"]   >= CUT) ? val["class"]   : ""
            ordr    = (("order"   in val) && boot["order"]   >= CUT) ? val["order"]   : ""
            family  = (("family"  in val) && boot["family"]  >= CUT) ? val["family"]  : ""
            genus   = (("genus"   in val) && boot["genus"]   >= CUT) ? val["genus"]   : ""
            species = (("species" in val) && boot["species"] >= CUT) ? val["species"] : ""
            gb = ("genus"   in boot) ? sprintf("%.4f", boot["genus"])   : ""
            sb = ("species" in boot) ? sprintf("%.4f", boot["species"]) : ""
            print seq, kingdom, phylum, cls, ordr, family, genus, species, gb, sb
        }
        ' ${prefix}_rdp_raw.txt > ${prefix}.taxonomy.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rdp_classifier: \$(rdp_classifier --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+' | head -1 || echo "2.13")
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
        // Handles SILVA, UNITE, MIDORI, coidb, PR2 formatted databases
        def add_species = (addspecies_db.name != 'NO_FILE') ? 'TRUE' : 'FALSE'

        """
        #!/usr/bin/env Rscript
        library(dada2)
        library(Biostrings)

        marker         <- "${marker}"
        db_type        <- "${db_type}"
        fasta_file     <- "${fasta}"
        tax_db         <- "${tax_db}"
        addspecies_db  <- "${addspecies_db}"
        prefix         <- "${prefix}"
        add_species    <- as.logical("${add_species}")

        # readDNAStringSet (unlike getSequences) preserves the FASTA headers
        # (ASV1, ASV2, ...) as names -- assignTaxonomy() below only accepts a
        # plain sequence vector and re-keys its own output by sequence,
        # discarding those names, so build a sequence->label map up front and
        # remap the result's rownames back afterward. This keeps asv_id as
        # the stable ASV label (matching the rdp/blast branches) instead of
        # the raw sequence, which downstream joins (asv_taxonomy_table) rely on.
        fasta_seqs   <- readDNAStringSet(fasta_file)
        seqs         <- as.character(fasta_seqs)
        seq_to_label <- setNames(names(seqs), unname(seqs))

        if (length(seqs) == 0) {
            message("Empty FASTA — writing empty taxonomy table.")
            write.table(data.frame(asv_id=character(0)), paste0(prefix, ".taxonomy.tsv"),
                        sep="\\t", quote=FALSE, row.names=FALSE, na="")
            writeLines(c(
                paste0('"${task.process}":'),
                paste0('    dada2: ', packageVersion('dada2')),
                paste0('    R: ', R.version\$major, '.', R.version\$minor)
            ), "versions.yml")
            quit(status=0)
        }

        # Assign taxonomy
        taxa <- assignTaxonomy(
            unname(seqs),
            tax_db,
            multithread  = ${task.cpus},
            minBoot      = 50,
            outputBootstraps = TRUE,
            tryRC        = TRUE
        )

        tax_table   <- taxa\$tax
        boot_table  <- taxa\$boot

        # addSpecies() does exact/near-exact sequence matching for species-level
        # ID -- much more reliable than naive-Bayes bootstrap at that depth.
        # Must run before the rowname remap below: it matches by the ASV's
        # actual sequence, which is what tax_table's rownames still are here.
        # No multithread param exists for addSpecies() -- it scans the
        # reference in chunks of `n` sequences regardless of query count, so
        # for large references (coidb is ~4.7M seqs) the chunk count (not
        # thread count) is what to tune; raise n well above the default 2000.
        if (add_species) {
            tax_table <- addSpecies(tax_table, addspecies_db, n=20000, verbose=TRUE)
        }

        rownames(tax_table)  <- unname(seq_to_label[rownames(tax_table)])
        rownames(boot_table) <- unname(seq_to_label[rownames(boot_table)])

        # Format ranks based on database type
        if (db_type %in% c("silva", "midori", "coidb")) {
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
    def prefix = (tax_method == 'dada2') ? "${marker}_${task.index}" : "${marker}"
    """
    touch ${prefix}.taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.3
    END_VERSIONS
    """

}
