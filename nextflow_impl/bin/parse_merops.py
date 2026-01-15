#!/usr/bin/env python3

import sys
import argparse
import logging
import csv
import re

def load_merops_mapping(lib_file):
    """
    Load pepunit.lib to map Sequence ID -> Family
    Header format: >mer_id #Family# ...
    """
    mapping = {}
    try:
        with open(lib_file, 'r', errors='replace') as f:
            for line in f:
                if line.startswith(">"):
                    # >MER00001 #Peptidase_S1# ...
                    # Remove \r if any (legacy tr/\015//d)
                    line = line.strip().replace('\r', '')
                    # Extract ID and Family
                    # Perl: ($mer_id) = $_ =~ /^>(.+?)\s/; -> Everything up to space
                    # But Diamond result uses truncated ID sometimes? Usually full ID.
                    # Diamond output format 6 (m8) usually preserves ID.
                    
                    parts = line.split(maxsplit=1)
                    if not parts: continue
                    full_id = parts[0][1:] # remove >
                    
                    # Search for #Family#
                    # Perl: $MEROPS_map{$tmp[1]} =~ /\#(.+?)\#/;
                    # So we need to store the header and extract family later?
                    # Or map ID -> Family now.
                    
                    match = re.search(r'\#(.+?)\#', line)
                    if match:
                        family = match.group(1)
                        mapping[full_id] = family
    except FileNotFoundError:
        pass # Should handle error 
    return mapping

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="m8 file")
    parser.add_argument("lib", help="pepunit.lib")
    parser.add_argument("output", help="output tsv")
    parser.add_argument("--genome_id", required=True)
    args = parser.parse_args()
    
    mapping = load_merops_mapping(args.lib)
    
    with open(args.input) as f, open(args.output, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["seq_id", "family", "target_id", "evalue", "genome_id"]) # Header
        
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split("\t")
            if len(parts) < 10: continue
            
            seq_id = parts[0]
            target_id = parts[1]
            evalue = parts[10] # m8: q, s, pident, len, mismatch, gap, qstart, qend, sstart, send, evalue, bitscore
            
            # Map target -> family
            family = mapping.get(target_id)
            if not family:
                # Try lookup without strict match?
                # Sometimes IDs diff?
                # Assume exact match for now.
                continue
                
            writer.writerow([seq_id, family, target_id, evalue, args.genome_id])

if __name__ == "__main__":
    main()
