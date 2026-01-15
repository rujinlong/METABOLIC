#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.raw_dir = "${projectDir}/tests/raw_db"
params.db_dir = "${projectDir}/tests/database"

process BUILD_HMM_KO {
    tag "kofam"
    publishDir "${params.db_dir}/kofam_database", mode: 'copy'

    input:
    path fasta

    output:
    path "kofam_merged.hmm"

    script:
    """
    hmmbuild --amino kofam_merged.hmm $fasta
    """
}

process BUILD_HMM_CUSTOM {
    tag "custom"
    publishDir "${params.db_dir}/METABOLIC_hmm_db", mode: 'copy'

    input:
    path fasta

    output:
    path "custom.hmm"

    script:
    """
    hmmbuild --amino custom.hmm $fasta
    """
}

process BUILD_HMM_DBCAN {
    tag "dbcan"
    publishDir "${params.db_dir}/dbCAN2", mode: 'copy'

    input:
    path fasta

    output:
    path "dbCAN-fam-HMMs.txt*"

    script:
    """
    hmmbuild --amino dbCAN-fam-HMMs.txt $fasta
    hmmpress dbCAN-fam-HMMs.txt
    """
}

process BUILD_DIAMOND {
    tag "merops"
    publishDir "${params.db_dir}/MEROPS", mode: 'copy'

    input:
    path fasta

    output:
    path "pepunit.db.dmnd"

    script:
    """
    diamond makedb --in $fasta --db pepunit.db
    """
}

workflow {
    ch_ko_fasta = Channel.fromPath("${params.raw_dir}/K00001.fasta")
    ch_cust_fasta = Channel.fromPath("${params.raw_dir}/custom.fasta")
    ch_dbcan_fasta = Channel.fromPath("${params.raw_dir}/GH1.fasta")
    ch_merops_fasta = Channel.fromPath("${params.raw_dir}/pepunit.lib")

    BUILD_HMM_KO(ch_ko_fasta)
    BUILD_HMM_CUSTOM(ch_cust_fasta)
    BUILD_HMM_DBCAN(ch_dbcan_fasta)
    BUILD_DIAMOND(ch_merops_fasta)
}
