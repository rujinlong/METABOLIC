process PARSE_HMM {
    tag "$genome_id"
    label 'process_low'
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    tuple val(genome_id), path(kofam_tbl), path(custom_tbl), path(proteins)
    path ko_list           // KOfam thresholds (K00001\tthreshold\tscore_type)
    path hmm_template      // Custom HMM thresholds (column 6=HMM, column 11=threshold)
    path motif_file
    path motif_pair_file

    output:
    tuple val(genome_id), path("${genome_id}_hits.tsv"), emit: hits

    script:
    def custom_arg = custom_tbl.name != 'NO_FILE' ? "--custom_results $custom_tbl" : ""
    def motif_arg = motif_file.name != 'NO_FILE' ? "--motif_file $motif_file" : ""
    def pair_arg = motif_pair_file.name != 'NO_FILE' ? "--motif_pair_file $motif_pair_file" : ""
    
    """
    python3 ${projectDir}/bin/parse_hmm_results.py \
        --kofam_results $kofam_tbl \
        $custom_arg \
        --proteins $proteins \
        --ko_list $ko_list \
        --hmm_template $hmm_template \
        $motif_arg \
        $pair_arg \
        --genome_id $genome_id \
        --output ${genome_id}_hits.tsv
    """
}

process PARSE_DBCAN {
    tag "$genome_id"
    label 'process_low'
    
    input:
    tuple val(genome_id), path(domtbl)

    output:
    path("${genome_id}_dbcan_parsed.tsv"), emit: hits

    script:
    """
    python3 ${projectDir}/bin/parse_dbcan.py \
        $domtbl \
        ${genome_id}_dbcan_parsed.tsv \
        --genome_id $genome_id
    """
}

process PARSE_MEROPS {
    tag "$genome_id"
    label 'process_low'
    
    input:
    tuple val(genome_id), path(m8)
    path merops_lib

    output:
    path("${genome_id}_merops_parsed.tsv"), emit: hits

    script:
    """
    python3 ${projectDir}/bin/parse_merops.py \
        $m8 \
        $merops_lib \
        ${genome_id}_merops_parsed.tsv \
        --genome_id $genome_id
    """
}
