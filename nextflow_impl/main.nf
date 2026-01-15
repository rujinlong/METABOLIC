#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Import Modules
 */
include { PRODIGAL } from './modules/annotation'
include { HMMSEARCH_KO; HMMSEARCH_CUSTOM } from './modules/search'
include { DBCAN_SEARCH } from './modules/dbcan'
include { MEROPS_SEARCH } from './modules/merops'
include { PARSE_HMM; PARSE_DBCAN; PARSE_MEROPS } from './modules/parsing'
include { GENERATE_TABLES } from './modules/tables'

/*
 * Help Message
 */
def helpMessage() {
    log.info """
    METABOLIC-Nextflow Pipeline
    ===========================
    Usage:
      nextflow run main.nf --input_genomes <path> [options]

    Required:
      --input_genomes     Path to directory containing genome FASTA files (*.fasta)

    Options:
      --outdir            Output directory (default: results)
      --kofam_dir         Path to merged KOfam HMM database
      --metabolic_hmm_dir Path to METABOLIC custom HMM directory
      --help              Show this help message

    Profiles:
      -profile standard   Run locally
      -profile docker     Run with Docker
      -profile singularity Run with Singularity
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.input_genomes) {
    log.error "Error: --input_genomes is required"
    helpMessage()
    exit 1
}

/*
 * Main Workflow
 */
workflow {
    
    // 1. Input Handling
    // Genomes: tuple(id, file)
    ch_genomes = Channel.fromPath("${params.input_genomes}/*.fasta")
                        .map { file -> tuple(file.simpleName, file) }

    // 2. Annotation
    PRODIGAL(ch_genomes)
    ch_proteins = PRODIGAL.out.proteins

    // 3. Search (HMM / dbCAN / MEROPS)
    
    // KOfam Search (Merged DB)
    ch_kofam_db = file(params.kofam_dir)
    
    // Custom HMM Search
    ch_custom_db = Channel.fromPath("${params.metabolic_hmm_dir}/*.hmm")
                          .filter { !(it.getFileName().toString() =~ /^K\d{5}\.hmm$/) }
                          .collectFile(name: 'custom_merged.hmm')
                          
    // dbCAN Search - need HMM + auxiliary h3* files
    ch_dbcan_db = Channel.fromPath("${params.dbcan_db}*").collect()
    
    // MEROPS Search
    ch_merops_db = file(params.merops_db)
    ch_merops_lib = file(params.merops_lib)
    
    // Run Searches
    HMMSEARCH_KO(ch_proteins, ch_kofam_db)
    HMMSEARCH_CUSTOM(ch_proteins, ch_custom_db)
    DBCAN_SEARCH(ch_proteins, ch_dbcan_db)
    MEROPS_SEARCH(ch_proteins, ch_merops_db)
    
    // Join HMM Results for Parsing
    ch_hmm_results = HMMSEARCH_KO.out.tblout
        .join(HMMSEARCH_CUSTOM.out.tblout, by: 0)
        .join(ch_proteins, by: 0)
        
    // 4. Parse Results
    
    // Parse HMM
    PARSE_HMM(
        ch_hmm_results,
        file(params.kofam_threshold),
        file(params.motif_file),
        file(params.motif_pair_file)
    )
    
    // Parse dbCAN
    PARSE_DBCAN(DBCAN_SEARCH.out.domtblout)
    
    // Parse MEROPS
    PARSE_MEROPS(MEROPS_SEARCH.out.m8, ch_merops_lib)
    
    // 5. Generate Worksheets
    // Collect all hits
    ch_all_hmm_hits = PARSE_HMM.out.hits.collectFile(name: 'all_hits_parsed.tsv', keepHeader: true)
    
    // Collect output from dbCAN/MEROPS
    ch_all_dbcan_hits = PARSE_DBCAN.out.hits.collectFile(name: 'all_dbcan_hits.tsv')
    ch_all_merops_hits = PARSE_MEROPS.out.hits.collectFile(name: 'all_merops_hits.tsv', keepHeader: true)
    
    // Genome IDs
    ch_genome_ids = ch_genomes.map { it[0] }.collect()
    
    GENERATE_TABLES(
        ch_all_hmm_hits,
        file(params.hmm_template_1),
        file(params.hmm_template_2),
        file(params.kegg_module_db),
        file(params.kegg_steps_db),
        ch_all_dbcan_hits,
        ch_all_merops_hits,
        ch_genome_ids
    )
}
