process SETUP_DATABASES {
    label 'process_single'
    executor 'local'

    input:
    val markers

    output:
    val true, emit: ready

    script:
    def marker_str = (markers instanceof List ? markers : [markers]).join(' ')
    """
    bash ${projectDir}/assets/download_databases.sh ${projectDir}/databases ${marker_str}
    """
}
