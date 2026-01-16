process PARSE_HMM {
    tag "all_proteins"
    label 'process_medium'
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    tuple val(id), path(kofam_tbl), path(custom_tbl), path(proteins)
    path ko_list
    path hmm_template
    path motif_file
    path motif_pair_file
    path contig_map
    path cluster_tsv

    output:
    path("all_hits_parsed.tsv"), emit: hits

    script:
    def custom_arg = custom_tbl.name != 'NO_FILE' ? "--custom_results $custom_tbl" : ""
    def motif_arg = motif_file.name != 'NO_FILE' ? "--motif_file $motif_file" : ""
    def pair_arg = motif_pair_file.name != 'NO_FILE' ? "--motif_pair_file $motif_pair_file" : ""
    def cluster_arg = cluster_tsv.name != 'NO_CLUSTER' ? "--cluster_tsv $cluster_tsv" : ""
    
    """
    python3 ${projectDir}/bin/parse_hmm_results.py \
        --kofam_results $kofam_tbl \
        $custom_arg \
        --proteins $proteins \
        --ko_list $ko_list \
        --hmm_template $hmm_template \
        $motif_arg \
        $pair_arg \
        --contig_map $contig_map \
        $cluster_arg \
        --output all_hits_parsed.tsv
    """
}

process PARSE_DBCAN {
    tag "all_proteins"
    label 'process_low'
    publishDir "${params.outdir}/tables", mode: 'copy'
    
    input:
    tuple val(id), path(domtbl)
    path contig_map
    path cluster_tsv

    output:
    path("all_dbcan_parsed.tsv"), emit: hits

    script:
    def cluster_arg = cluster_tsv.name != 'NO_CLUSTER' ? "--cluster_tsv $cluster_tsv" : ""
    
    """
    python3 ${projectDir}/bin/parse_dbcan.py \
        $domtbl \
        all_dbcan_parsed.tsv \
        --contig_map $contig_map \
        $cluster_arg
    """
}

process PARSE_MEROPS {
    tag "all_proteins"
    label 'process_low'
    publishDir "${params.outdir}/tables", mode: 'copy'
    
    input:
    tuple val(id), path(m8)
    path merops_lib
    path contig_map
    path cluster_tsv

    output:
    path("all_merops_parsed.tsv"), emit: hits

    script:
    def cluster_arg = cluster_tsv.name != 'NO_CLUSTER' ? "--cluster_tsv $cluster_tsv" : ""
    
    """
    python3 ${projectDir}/bin/parse_merops.py \
        $m8 \
        $merops_lib \
        all_merops_parsed.tsv \
        --contig_map $contig_map \
        $cluster_arg
    """
}
