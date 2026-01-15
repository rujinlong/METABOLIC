#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import re
from Bio import SeqIO

def load_thresholds(threshold_file, db_type='full'):
    """
    Load KOfam thresholds. 
    Format: K00001 \t threshold \t score_type
    Returns dict: {'K00001.hmm': {'threshold': 100.0, 'type': 'full'}}
    """
    thresholds = {}
    with open(threshold_file) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                # Legacy code handles '-' as a specific case, usually default
                hmm = f"{parts[0]}.hmm"
                if parts[1] == "-":
                     thresholds[hmm] = {'threshold': 50.0, 'type': 'full'}
                else:
                    thresholds[hmm] = {'threshold': float(parts[1]), 'type': parts[2]}
    return thresholds

def load_motifs(motif_file):
    """
    Load motif regex definitions.
    Format: hmm_name:motif_string (with X)
    Returns dict: {'hmm': compiled_regex}
    """
    motifs = {}
    if not motif_file: return motifs
    
    with open(motif_file) as f:
        for line in f:
            parts = line.strip().split(":")
            if len(parts) == 2:
                hmm, pattern = parts
                # Convert X to amino acid character set
                pattern = pattern.replace("X", "[ARNDCQEGHILKMFPSTWYV]")
                motifs[hmm] = re.compile(pattern)
    return motifs

def load_motif_pairs(pair_file):
    """
    Load motif pairs. Format: hmm:partner_hmm
    """
    pairs = {}
    if not pair_file: return pairs
    with open(pair_file) as f:
        for line in f:
            parts = line.strip().split(":")
            if len(parts) == 2:
                pairs[parts[0]] = parts[1]
    return pairs

def parse_tblout(tblout_file):
    """
    Parse generic tblout file into DataFrame.
    Columns: seq_id, accession, hmm_name, accession, ..., full_score, ..., domain_score
    """
    # tblout is whitespace delimited, fix widths? No, usually generic split is safer if no internal spaces
    # But usually extract: col 0 (target), col 2 (query/hmm), col 5 (full score), col 8 (best dom score)
    data = []
    with open(tblout_file) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split()
            if len(parts) < 10: continue
            data.append({
                'seq_id': parts[0],
                'hmm_name': parts[2],
                'full_score': float(parts[5]),
                'domain_score': float(parts[8])
            })
    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description="Parse HMMsearch results for METABOLIC")
    parser.add_argument("--kofam_results", required=True)
    parser.add_argument("--custom_results", required=False) # Can be merged or separate
    parser.add_argument("--proteins", required=True, help="Fasta file for motif validation")
    parser.add_argument("--kofam_thresholds", required=True)
    parser.add_argument("--motif_file", required=False)
    parser.add_argument("--motif_pair_file", required=False)
    parser.add_argument("--genome_id", required=True, help="Genome ID to attach to output")
    parser.add_argument("--output", required=True)
    
    args = parser.parse_args()

    # 1. Load Definitions
    thresholds = load_thresholds(args.kofam_thresholds)
    motif_regex = load_motifs(args.motif_file) if args.motif_file else {}
    motif_pairs = load_motif_pairs(args.motif_pair_file) if args.motif_pair_file else {}
    
    # 2. Parse Results (KOfam)
    df_ko = parse_tblout(args.kofam_results)
    
    # Merge custom results if exists
    if args.custom_results:
        df_cust = parse_tblout(args.custom_results)
        df = pd.concat([df_ko, df_cust])
    else:
        df = df_ko

    if df.empty:
        open(args.output, 'w').close()
        return

    # 3. Filter by Thresholds
    valid_hits = []
    
    # Pre-load sequences ONLY if motifs are required (Optimisation)
    # Check if any candidate hits map to HMMs that have motifs
    # Actually, legacy loads all. We can load lazy or load all.
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.proteins, "fasta"))

    for _, row in df.iterrows():
        hmm = f"{row['hmm_name']}.hmm"
        
        # Check thresholds
        if hmm not in thresholds:
            # Maybe skip or keep? Legacy says "next unless exists"
            continue
            
        th = thresholds[hmm]
        pass_threshold = False
        
        if th['type'] == 'domain':
            if row['domain_score'] >= th['threshold']:
                pass_threshold = True
        else:
            if row['full_score'] >= th['threshold']:
                pass_threshold = True
                
        if not pass_threshold:
            continue
            
        # 4. Filter by Motif / Pair
        row_id = (row['seq_id'], row['hmm_name'])
        
        if row['hmm_name'] in motif_regex:
            # Check sequence
            seq_record = seq_dict.get(row['seq_id'])
            if seq_record:
                seq_str = str(seq_record.seq)
                if motif_regex[row['hmm_name']].search(seq_str):
                     valid_hits.append(row)
        
        elif row['hmm_name'] in motif_pairs:
            # Pair logic: Postpone to batch processing
            # We need to collect all hits for this seq_id to compare scores later
            pass 
        else:
            # Plain hit for HMMs that are NOT in motif_regex AND NOT in motif_pairs
            # Wait, what if it IS a partner in a pair?
            # The motif_pairs dict is {hmm: partner}.
            # We need to know if 'hmm_name' acts as a partner too? 
            # In legacy: "if exists Motif_pair{$hmm}".
            # So if it's the KEY in the pair map, we check.
            # If it's the VALUE, what happens? 
            # Legacy only checks `if (exists $Motif_pair{$hmm_basename})`.
            # If it's the 'anti' hmm, it proceeds to 'else' block -> _process_hit.
            # BUT, the `push @motif_pair_candidates` logic eventually runs hmmsearch for BOTH.
            # And then: `if ($motif_score >= $anti_score ...)`
            # So effectively, for the PRIMARY motif (key), we check against secondary.
            # For the SECONDARY motif (value), does it have its own entry?
            # Usually pairs are directional A->B.
            valid_hits.append(row)

    # 5. Handle Motif Pairs
    # We need to handle the case where we saw the Key but need to check the Value's score.
    # The dataframe `df` contains ALL hits.
    # We need a lookup for scores: (seq_id, hmm_name) -> full_score
    
    # Create score lookup
    score_map = {} # (seq_id, hmm_name) -> score
    for _, row in df.iterrows():
         score_map[(row['seq_id'], row['hmm_name'])] = row['full_score']
    
    # Iterate thresholds again for pairs?
    # Or just iterate the original df and filter?
    # We only care about rows where `hmm_name` is in `motif_pairs`.
    
    for _, row in df.iterrows():
        hmm = row['hmm_name']
        if hmm in motif_pairs:
            # This is a Motif Pair Key
            partner_hmm = motif_pairs[hmm]
            
            my_score = row['full_score']
            
            # Check partner score
            partner_score = score_map.get((row['seq_id'], partner_hmm), 0.0)
            
            # Check threshold (already done implicitly? No, we skipped appending validity above)
            # Threshold check for 'hmm'
            full_hmm = f"{hmm}.hmm"
            if full_hmm in thresholds:
                 th = thresholds[full_hmm]
                 pass_th = False
                 if th['type'] == 'domain':
                      if row['domain_score'] >= th['threshold']: pass_th = True
                 else:
                      if row['full_score'] >= th['threshold']: pass_th = True
                 
                 if pass_th:
                     # Compare with partner
                     if my_score >= partner_score:
                         valid_hits.append(row)

    final_hits = pd.DataFrame(valid_hits)
    
    # Write Output
    if not final_hits.empty:
        # Deduplicate? If needed.
        # Legacy might output multiple hits.
        final_hits['genome_id'] = args.genome_id
        final_hits.to_csv(args.output, sep="\t", index=False)
    else:
        open(args.output, 'w').close()

if __name__ == "__main__":
    main()
