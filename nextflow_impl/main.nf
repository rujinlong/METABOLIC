#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Import Modules
 */
include { MMSEQS_CLUSTER } from './modules/cluster'
include { HMMSEARCH_KO; HMMSEARCH_CUSTOM } from './modules/search'
include { DBCAN_SEARCH } from './modules/dbcan'
include { MEROPS_SEARCH } from './modules/merops'
include { PARSE_HMM; PARSE_DBCAN; PARSE_MEROPS } from './modules/parsing'
include { GENERATE_TABLES } from './modules/tables'
include { GENERATE_R_INPUTS; PLOT_CYCLES; CREATE_EXCEL } from './modules/visualization'

/*
 * Help Message
 */
def helpMessage() {
    log.info """
    METABOLIC-Nextflow Pipeline
    ===========================
    
    A streamlined, HPC-optimized pipeline for metabolic potential profiling.

    Usage:
      nextflow run main.nf --proteins <fasta> --contig_map <tsv> [options]

    Required Inputs:
      --proteins        Path to a SINGLE concatenated protein FASTA file
                        Header format MUST be: >contigID_proteinID [description]
                        Example: >contig001_1 hypothetical protein
      
      --contig_map      Path to TSV file mapping contig IDs to genome/MAG IDs
                        Format: contigID<TAB>genomeID (no header)

    Clustering Options:
      --skip_cluster    Skip MMseqs2 clustering (default: false)
                        Use if input is already dereplicated.
      --cluster_sid     Minimum sequence identity for clustering (default: 0.9)
      --cluster_cov     Minimum coverage for clustering (default: 0.8)

    Other Options:
      --outdir          Output directory (default: results)
      --db_dir          Root directory for all METABOLIC databases
      --help            Show this help message

    Profiles:
      -profile standard     Run locally
      -profile docker       Run with Docker
      -profile singularity  Run with Singularity
      -profile hpc          Run on HPC with SLURM
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.proteins) {
    log.error "Error: --proteins is required"
    helpMessage()
    exit 1
}

if (!params.contig_map) {
    log.error "Error: --contig_map is required"
    helpMessage()
    exit 1
}

/*
 * Main Workflow
 */
workflow {
    
    // 1. Input: Single concatenated protein file
    log.info "Using protein sequences from: ${params.proteins}"
    log.info "Using contig-to-genome mapping: ${params.contig_map}"
    
    ch_input_proteins = file(params.proteins)
    ch_contig_map = file(params.contig_map)

    // 2. Optional Clustering Step
    if (params.skip_cluster) {
        log.info "Skipping MMseqs2 clustering (--skip_cluster=true)"
        ch_proteins = Channel.of(tuple('all_proteins', ch_input_proteins))
        ch_cluster_tsv = Channel.of(file('NO_CLUSTER'))
    } else {
        log.info "Running MMseqs2 clustering (identity=${params.cluster_sid}, coverage=${params.cluster_cov})"
        MMSEQS_CLUSTER(ch_input_proteins)
        ch_proteins = MMSEQS_CLUSTER.out.rep_seq.map { tuple('rep_proteins', it) }
        ch_cluster_tsv = MMSEQS_CLUSTER.out.cluster_tsv
    }

    // 3. Database channels
    ch_kofam_db = file(params.kofam_dir)
    ch_custom_db = file(params.metabolic_hmm)
    ch_dbcan_db = Channel.fromPath("${params.dbcan_db}*").collect()
    ch_merops_db = file(params.merops_db)
    ch_merops_lib = file(params.merops_lib)
    
    // 4. Run Searches (all in parallel)
    HMMSEARCH_KO(ch_proteins, ch_kofam_db)
    HMMSEARCH_CUSTOM(ch_proteins, ch_custom_db)
    DBCAN_SEARCH(ch_proteins, ch_dbcan_db)
    MEROPS_SEARCH(ch_proteins, ch_merops_db)
    
    // 5. Join HMM Results for Parsing
    ch_hmm_results = HMMSEARCH_KO.out.tblout
        .join(HMMSEARCH_CUSTOM.out.tblout, by: 0)
        .combine(ch_proteins.map { it[1] })  // Add protein file for motif validation
        
    // 6. Parse Results with contig_map and cluster_tsv for expansion
    PARSE_HMM(
        ch_hmm_results,
        file(params.ko_list),
        file(params.hmm_template_1),
        file(params.motif_file),
        file(params.motif_pair_file),
        ch_contig_map,
        ch_cluster_tsv
    )
    
    PARSE_DBCAN(DBCAN_SEARCH.out.domtblout, ch_contig_map, ch_cluster_tsv)
    
    PARSE_MEROPS(MEROPS_SEARCH.out.m8, ch_merops_lib, ch_contig_map, ch_cluster_tsv)
    
    // 7. Generate Worksheets
    ch_genome_ids = Channel.fromPath(params.contig_map)
        .splitCsv(sep: '\t')
        .map { it[1] }
        .unique()
        .collect()
    
    GENERATE_TABLES(
        PARSE_HMM.out.hits,
        file(params.hmm_template_1),
        file(params.hmm_template_2),
        file(params.kegg_module_db),
        file(params.kegg_steps_db),
        PARSE_DBCAN.out.hits,
        PARSE_MEROPS.out.hits,
        ch_genome_ids
    )
    
    // 8. Generate Visualization (R_input files + PDF cycle diagrams)
    GENERATE_R_INPUTS(
        PARSE_HMM.out.hits,
        file(params.r_pathways)
    )
    
    PLOT_CYCLES(GENERATE_R_INPUTS.out.r_input_dir)
    
    // 9. Create Excel spreadsheet from TSV worksheets
    CREATE_EXCEL(GENERATE_TABLES.out.tables_dir)
}
