process SINTAX_BLANK_CLEANUP {
    tag "${marker}"
    label 'process_single'

    container 'ghcr.io/ukceh-molecularecology/nf-edna-pipeline/tapirs-python:1.0'

    publishDir "${params.outdir}/sintax/${marker}/blank_cleanup", mode: 'copy'

    input:
    tuple val(marker), path(clare_abundance_csv), path(ncl_matrix_csv)

    output:
    tuple val(marker), path('ncl_cleaned_bothLOD.csv'), emit: cleaned_wide
    tuple val(marker), path('ncl_cleaned_labLOD.csv'),  emit: cleaned_lab
    tuple val(marker), path('ncl_cleaned_siteLOD.csv'), emit: cleaned_site
    tuple val(marker), path('ncl_cleaned_long.csv'),    emit: cleaned_long
    tuple val(marker), path('cleanup_*.csv'),           emit: diagnostics
    tuple val(marker), path('*_pa.csv'),                emit: presence_absence
    path 'versions.yml',                                emit: versions

    script:
    def bc = params.sintax.blank_cleanup
    def exclusion_flag = bc.enable_taxon_exclusion ? '--enable-taxon-exclusion' : ''
    def excluded_taxa   = bc.excluded_taxa.join(',')
    """
    blank_cleanup.py \\
        --clare-abundance ${clare_abundance_csv} \\
        --ncl-matrix ${ncl_matrix_csv} \\
        --pcr-blank-regex '${bc.pcr_blank_regex}' \\
        --extraction-blank-regex '${bc.extraction_blank_regex}' \\
        --site-blank-regex '${bc.site_blank_regex}' \\
        --lod-sd-multiplier ${bc.lod_sd_multiplier} \\
        --loq-sd-multiplier ${bc.loq_sd_multiplier} \\
        --excluded-taxa '${excluded_taxa}' \\
        ${exclusion_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ncl_cleaned_bothLOD.csv ncl_cleaned_labLOD.csv ncl_cleaned_siteLOD.csv ncl_cleaned_long.csv
    touch cleanup_summary.csv
    touch ncl_cleaned_bothLOD_pa.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
