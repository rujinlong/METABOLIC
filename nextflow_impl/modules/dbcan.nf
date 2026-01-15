process DBCAN_SEARCH {
    tag "$id"
    label 'process_high'
    // publishDir "${params.outdir}/dbcan", mode: 'copy'

    input:
    tuple val(id), path(proteins)
    path dbcan_db

    output:
    tuple val(id), path("${id}.dbcan.dm"), emit: domtblout

    script:
    """
    hmmscan --domtblout ${id}.dbcan.dm \
            --cpu ${task.cpus} \
            $dbcan_db \
            $proteins
    """
}
