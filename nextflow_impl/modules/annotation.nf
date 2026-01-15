process PRODIGAL {
    tag "$genome_id"
    label 'process_medium'
    publishDir "${params.outdir}/proteins", mode: 'copy'

    input:
    tuple val(genome_id), path(fasta)

    output:
    tuple val(genome_id), path("${genome_id}.faa"), emit: proteins
    tuple val(genome_id), path("${genome_id}.gff"), emit: gff

    script:
    """
    prodigal -i $fasta \
             -a ${genome_id}.faa \
             -o ${genome_id}.gff \
             -p meta -f gff -q
    """
}
