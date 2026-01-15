process DBCAN_SEARCH {
    tag "$id"
    label 'process_high'
    // publishDir "${params.outdir}/dbcan", mode: 'copy'

    input:
    tuple val(id), path(proteins)
    path dbcan_files  // All files: HMM + .h3f, .h3i, .h3m, .h3p

    output:
    tuple val(id), path("${id}.dbcan.dm"), emit: domtblout

    script:
    // Find the main HMM file (the one without .h3 extension)
    def hmm_file = dbcan_files.find { it.name.endsWith('.txt') || it.name.endsWith('.hmm') }
    """
    hmmscan --domtblout ${id}.dbcan.dm \
            --cpu ${task.cpus} \
            ${hmm_file} \
            $proteins
    """
}
