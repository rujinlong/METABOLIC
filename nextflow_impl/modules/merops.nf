process MEROPS_SEARCH {
    tag "$id"
    label 'process_medium'
    // publishDir "${params.outdir}/merops", mode: 'copy'

    input:
    tuple val(id), path(proteins)
    path merops_db

    output:
    tuple val(id), path("${id}.merops.m8"), emit: m8

    script:
    // Diamond expects DB name without .dmnd extension
    def db_name = merops_db.baseName.replaceAll(/\.dmnd$/, '')
    """
    diamond blastp \
        -d $db_name \
        -q $proteins \
        -o ${id}.merops.m8 \
        -k 1 \
        -e 1e-10 \
        --query-cover 80 \
        --id 50 \
        --quiet \
        -p ${task.cpus}
    """
}
