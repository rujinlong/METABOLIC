#!/usr/bin/env python3

import os
import sys

def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def write_file(path, content):
    with open(path, 'w') as f:
        f.write(content)

def main():
    root = "tests"
    create_dir(root)
    
    # 1. Input Genomes
    genome_dir = os.path.join(root, "data")
    create_dir(genome_dir)
    
    # DNA sequence for Prodigal (Nucleotides)
    # Met (ATG) - Ala (GCA) * 60 - Stop (TAA)
    # Prodigal needs DNA input to predict proteins
    gene1 = "ATG" + "GCA"*60 + "TAA"
    gene2 = "ATG" + "GCC"*80 + "TGA"
    
    dna_seq = gene1 + "NNNNNNNNNN" + gene2
    
    # Write Genome (DNA)
    genome_fasta = f">contig1\n{dna_seq}\n"
    write_file(os.path.join(genome_dir, "test_genome.fasta"), genome_fasta)
    
    # Protein sequence for HMM building
    protein_seq = "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    
    # 2. Database Build Inputs (Raw FASTAs for hmmbuild/diamond)
    raw_db_dir = os.path.join(root, "raw_db")
    create_dir(raw_db_dir)
    
    # Alignment for K00001 (Needs multiple seqs for hmmbuild usually, or 1 works)
    # 3 sequences
    align_fasta = f">s1\n{protein_seq}\n>s2\n{protein_seq}\n>s3\n{protein_seq}\n"
    write_file(os.path.join(raw_db_dir, "K00001.fasta"), align_fasta)
    
    # Alignment for Custom
    write_file(os.path.join(raw_db_dir, "custom.fasta"), align_fasta)
    
    # Alignment for dbCAN (GH1)
    write_file(os.path.join(raw_db_dir, "GH1.fasta"), align_fasta)
    
    # FASTA for MEROPS (pepunit.lib)
    # Needs special header for parsing: >ID #Family#
    merops_fasta = f">MER00001 #Peptidase_A1# Test peptidase\n{protein_seq}\n"
    write_file(os.path.join(raw_db_dir, "pepunit.lib"), merops_fasta)
    
    # 3. Template/Ancillary Files
    # Note: These need to be in the final 'database' folder structure
    db_root = os.path.join(root, "database")
    create_dir(db_root)
    
    # METABOLIC_template_and_database
    meta_tmpl_dir = os.path.join(db_root, "METABOLIC_template_and_database")
    create_dir(meta_tmpl_dir)
    
    # prokaryote.hal
    # K00001 \t threshold \t type
    # Using low threshold to ensure hit
    write_file(os.path.join(meta_tmpl_dir, "prokaryote.hal"), "K00001\t10.0\tfull\n")
    
    # motif.txt (Empty for now or dummy)
    write_file(os.path.join(meta_tmpl_dir, "motif.txt"), "")
    
    # motif.pair.txt
    write_file(os.path.join(meta_tmpl_dir, "motif.pair.txt"), "")
    
    # hmm_table_template.txt
    # Entry \t Category \t Function \t GeneAbbr \t GeneName \t HmmFile \t KO \t Reaction \t ... \t Product \t Threshold
    # Map K00001
    tmpl1 = "Entry1\tMetabolism\tTestFunc\tgeneA\tGeneA\tK00001.hmm,custom.hmm\tK00001\tR00001\t\tProductA\t10.0\n"
    write_file(os.path.join(meta_tmpl_dir, "hmm_table_template.txt"), tmpl1)
    
    # hmm_table_template_2.txt
    # ID \t EntryRefs \t Cat \t Func \t Abbr
    tmpl2 = "Func1\tEntry1\tMetabolism\tTest Function Group\tF1\n"
    write_file(os.path.join(meta_tmpl_dir, "hmm_table_template_2.txt"), tmpl2)
    
    # ko00002.keg
    # C    Category
    # D      M00001  Desc
    keg_content = "C    Metabolism\nD      M00001  Test Module\n"
    write_file(os.path.join(meta_tmpl_dir, "ko00002.keg"), keg_content)
    
    # kegg_module_step_db.txt
    # name \t K-string \t StepID
    # Desc \t (K00001) \t M00001+01
    steps_content = "name\tdef\tstep\nTest Step\tK00001\tM00001+01\n"
    write_file(os.path.join(meta_tmpl_dir, "kegg_module_step_db.txt"), steps_content)
    
    # 4. Create placeholders for DB dirs 
    # (setup_test.nf will fill them)
    create_dir(os.path.join(db_root, "kofam_database"))
    create_dir(os.path.join(db_root, "METABOLIC_hmm_db"))
    create_dir(os.path.join(db_root, "dbCAN2"))
    create_dir(os.path.join(db_root, "MEROPS"))
    
    # Copy pepunit.lib to MEROPS dir as it is needed as text too
    write_file(os.path.join(db_root, "MEROPS", "pepunit.lib"), merops_fasta)

    print("Created test data inputs in tests/")
    print("Run 'nextflow run setup_test.nf' to build binaries.")

if __name__ == "__main__":
    main()
