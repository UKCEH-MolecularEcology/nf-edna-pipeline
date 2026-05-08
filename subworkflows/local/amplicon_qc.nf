include { FASTQC   } from '../../modules/local/fastqc/main'
include { CUTADAPT } from '../../modules/local/cutadapt/main'

workflow AMPLICON_QC {

    take:
    ch_reads    // [ meta, reads ]

    main:
    ch_versions = Channel.empty()

    // Raw read QC
    FASTQC(ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Primer trimming (marker-specific primers embedded in meta from nextflow.config)
    CUTADAPT(ch_reads)
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    emit:
    trimmed_reads = CUTADAPT.out.reads
    fastqc_zip    = FASTQC.out.zip
    cutadapt_log  = CUTADAPT.out.log
    versions      = ch_versions
}
