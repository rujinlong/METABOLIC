process PYRODIGAL {
    tag "$genome_id"
    label 'process_high'
    publishDir "${params.outdir}/proteins", mode: 'copy'

    input:
    tuple val(genome_id), path(fasta)

    output:
    tuple val(genome_id), path("${genome_id}.faa"), emit: proteins
    tuple val(genome_id), path("${genome_id}.gff"), emit: gff
    tuple val(genome_id), path("${genome_id}.fna"), emit: nucleotides

    script:
    """
    pyrodigal -i $fasta \
              -a ${genome_id}.faa \
              -d ${genome_id}.fna \
              -o ${genome_id}.gff \
              -p ${params.pyrodigal_mode} \
              -f gff \
              -m \
              --pool process \
              -j ${task.cpus}
    """
}
