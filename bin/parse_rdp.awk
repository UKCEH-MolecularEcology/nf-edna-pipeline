BEGIN {
    FS = "\t"
    prefix["domain"]       = "k__"
    prefix["superkingdom"] = "k__"
    prefix["kingdom"]      = "k__"
    prefix["phylum"]       = "p__"
    prefix["class"]        = "c__"
    prefix["order"]        = "o__"
    prefix["family"]       = "f__"
    prefix["genus"]        = "g__"
    prefix["species"]      = "s__"
    print "ASV_ID\ttaxonomy\tconfidence"
    n_total = 0; n_classified = 0
}
NF == 0 { next }
substr($0, 1, 1) == "#" { next }
NF < 3 { next }
{
    seq_id = $1
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", seq_id)
    if (seq_id == "") next
    n_total++
    # Field 2 is always the orientation field ('+', '-', or empty).
    # Triplets (name, rank, conf) always start at field 3.
    start = 3
    taxonomy = ""; last_conf = 1.0
    for (i = start; i + 2 <= NF; i += 3) {
        name = $i; rank = $(i+1); conf = $(i+2) + 0
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", name)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", rank)
        rank = tolower(rank)
        if (rank == "rootrank" || name == "" || tolower(name) == "root") continue
        if (!(rank in prefix)) continue
        if (conf < threshold) break
        taxonomy = (taxonomy == "") ? prefix[rank] name : taxonomy ";" prefix[rank] name
        last_conf = conf
    }
    if (taxonomy == "") taxonomy = "k__unclassified"
    else n_classified++
    printf "%s\t%s\t%.4f\n", seq_id, taxonomy, last_conf
}
END {
    print "Taxonomy: " n_classified "/" n_total " sequences classified above threshold " threshold > "/dev/stderr"
}
