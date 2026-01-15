process HMMSEARCH {
    tag "$id"
    label 'process_high'
    // Don't publish intermediate tables unless debug
    // publishDir "${params.outdir}/hmm_results", mode: 'copy'

    input:
    tuple val(id), path(proteins)
    path hmm_db

    output:
    tuple val(id), path("${id}.tblout"), emit: tblout
    tuple val(id), path("${id}.domtblout"), emit: domtblout

    script:
    """
    hmmsearch --cpu ${task.cpus} \
              --noali \
              --tblout ${id}.tblout \
              --domtblout ${id}.domtblout \
              $hmm_db \
              $proteins
    """
}
