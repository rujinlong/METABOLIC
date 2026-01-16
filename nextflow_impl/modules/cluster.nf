process MMSEQS_CLUSTER {
    tag "clustering"
    label 'process_high'
    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path proteins

    output:
    path "rep_seq.fasta", emit: rep_seq
    path "cluster.tsv", emit: cluster_tsv

    script:
    """
    # Create MMseqs2 database
    mmseqs createdb $proteins proteins_db
    
    # Cluster proteins
    mmseqs cluster proteins_db cluster_db tmp \
        --min-seq-id ${params.cluster_sid} \
        -c ${params.cluster_cov} \
        --cov-mode 0 \
        --threads ${task.cpus}
    
    # Extract representative sequences
    mmseqs createsubdb cluster_db proteins_db rep_db
    mmseqs convert2fasta rep_db rep_seq.fasta
    
    # Create cluster membership TSV: representative_id <TAB> member_id
    mmseqs createtsv proteins_db proteins_db cluster_db cluster.tsv
    
    # Clean up
    rm -rf tmp
    """
}
