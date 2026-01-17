// Process to split FASTA into chunks for parallel processing

process SPLIT_FASTA {
    tag "split"
    label 'process_low'
    
    input:
    tuple val(name), path(fasta)
    val chunk_size
    
    output:
    path "chunk_*.faa", emit: chunks
    
    script:
    """
    split_fasta.py ${fasta} --chunk_size ${chunk_size} --prefix chunk_
    """
}
