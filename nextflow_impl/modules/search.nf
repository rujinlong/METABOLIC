process HMMSEARCH_KO {
    tag "$id"
    label 'process_high'

    input:
    tuple val(id), path(proteins)
    path hmm_db

    output:
    tuple val(id), path("${id}.kofam.tblout"), emit: tblout
    tuple val(id), path("${id}.kofam.domtblout"), emit: domtblout

    script:
    """
    hmmsearch --cpu ${task.cpus} \
              --noali \
              --tblout ${id}.kofam.tblout \
              --domtblout ${id}.kofam.domtblout \
              $hmm_db \
              $proteins
    """
}

process HMMSEARCH_CUSTOM {
    tag "$id"
    label 'process_high'

    input:
    tuple val(id), path(proteins)
    path hmm_db

    output:
    tuple val(id), path("${id}.custom.tblout"), emit: tblout
    tuple val(id), path("${id}.custom.domtblout"), emit: domtblout

    script:
    """
    hmmsearch --cpu ${task.cpus} \
              --noali \
              --tblout ${id}.custom.tblout \
              --domtblout ${id}.custom.domtblout \
              $hmm_db \
              $proteins
    """
}
