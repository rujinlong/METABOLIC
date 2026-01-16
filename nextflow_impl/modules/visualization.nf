// Visualization and Excel output processes

process GENERATE_R_INPUTS {
    tag "r_inputs"
    label 'process_low'
    publishDir "${params.outdir}/figures_input", mode: 'copy'

    input:
    path hits_tsv
    path r_pathways

    output:
    path "Nutrient_Cycling_Diagram_Input", emit: r_input_dir

    script:
    """
    generate_r_inputs.py ${hits_tsv} ${r_pathways} --outdir Nutrient_Cycling_Diagram_Input
    """
}


process PLOT_CYCLES {
    tag "cycles"
    label 'process_low'
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path r_input_dir

    output:
    path "draw_biogeochem_cycles/*.pdf", optional: true

    script:
    """
    Rscript ${params.metabolic_dir}/draw_biogeochemical_cycles.R \
        ${r_input_dir} Output FALSE > /dev/null 2>&1 || true
    
    # Move output if created
    if [ -d "Output/draw_biogeochem_cycles" ]; then
        mv Output/draw_biogeochem_cycles .
    else
        mkdir -p draw_biogeochem_cycles
        echo "No plots generated" > draw_biogeochem_cycles/README.txt
    fi
    """
}


process CREATE_EXCEL {
    tag "excel"
    label 'process_low'
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path tables_dir

    output:
    path "METABOLIC_result.xlsx", optional: true

    script:
    """
    # Copy TSV files to current directory for R script
    cp ${tables_dir}/METABOLIC_result_worksheet*.tsv .
    
    Rscript ${params.metabolic_dir}/create_excel_spreadsheet.R ./ > /dev/null 2>&1 || true
    
    # Check if Excel was created
    if [ ! -f "METABOLIC_result.xlsx" ]; then
        echo "Excel generation failed - R package openxlsx may be missing" >&2
    fi
    """
}
