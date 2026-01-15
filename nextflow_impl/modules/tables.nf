process GENERATE_TABLES {
    label 'process_low'
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path hits_file
    path template1
    path template2
    path kegg_module_db
    path kegg_steps
    path dbcan_hits
    path merops_hits
    val genome_ids

    output:
    path "METABOLIC_result_worksheet1.tsv"
    path "METABOLIC_result_worksheet2.tsv"
    path "METABOLIC_result_worksheet3.tsv"
    path "METABOLIC_result_worksheet4.tsv"
    path "METABOLIC_result_worksheet5.tsv"
    path "METABOLIC_result_worksheet6.tsv"

    script:
    def dbcan_arg = dbcan_hits.name != 'NO_FILE' ? "--dbcan_hits $dbcan_hits" : ""
    def merops_arg = merops_hits.name != 'NO_FILE' ? "--merops_hits $merops_hits" : ""

    """
    python3 ${projectDir}/bin/generate_worksheets.py \
        --hits $hits_file \
        --template1 $template1 \
        --template2 $template2 \
        --kegg_module_db $kegg_module_db \
        --kegg_steps $kegg_steps \
        $dbcan_arg \
        $merops_arg \
        --output_dir ./ \
        --genome_ids ${genome_ids.join(',')}
    """
}
