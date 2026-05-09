include { FASTQC } from '../../modules/local/fastqc/main'

workflow AMPLICON_QC {

    take:
    ch_reads    // [ meta, reads ]

    main:
    ch_versions = Channel.empty()

    // Raw read QC
    FASTQC(ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    emit:
    fastqc_zip = FASTQC.out.zip
    versions   = ch_versions
}
