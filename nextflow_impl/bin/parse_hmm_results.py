#!/usr/bin/env python3
"""
Parse HMMsearch results for METABOLIC.

This script filters HMM hits based on:
1. Per-HMM thresholds (from ko_list and hmm_table_template)
2. Motif validation (active site residues)
3. Motif pair comparisons (competing HMMs)

It uses a contig_map to resolve protein IDs to genome/MAG IDs.
"""

import sys
import re
from pathlib import Path
from typing import Optional

import polars as pl
import typer
from Bio import SeqIO
from rich.console import Console
from rich.progress import Progress

console = Console(stderr=True)
app = typer.Typer(add_completion=False)


def load_contig_map(contig_map_file: Path) -> dict[str, str]:
    """Load contig-to-genome mapping."""
    contig_to_genome = {}
    with open(contig_map_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                contig_to_genome[parts[0]] = parts[1]
    return contig_to_genome


def extract_contig_from_protein(protein_id: str) -> str:
    """Extract contigID from proteinID (format: contigID_proteinNumber)."""
    parts = protein_id.rsplit('_', 1)
    return parts[0] if len(parts) > 1 else protein_id


def load_cluster_map(cluster_tsv_file: Path) -> dict[str, list[str]]:
    """Load MMseqs2 cluster membership: rep_id -> [member_ids]."""
    cluster_map: dict[str, list[str]] = {}
    with open(cluster_tsv_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                rep_id, member_id = parts[0], parts[1]
                if rep_id not in cluster_map:
                    cluster_map[rep_id] = []
                cluster_map[rep_id].append(member_id)
    return cluster_map


def load_ko_list_thresholds(ko_list_file: Path) -> dict[str, dict]:
    """Load KOfam thresholds from ko_list file."""
    thresholds = {}
    with open(ko_list_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0].startswith('K'):
                hmm = f"{parts[0]}.hmm"
                if parts[1] == "-":
                    thresholds[hmm] = {'threshold': 50.0, 'type': 'full'}
                else:
                    try:
                        thresholds[hmm] = {'threshold': float(parts[1]), 'type': parts[2]}
                    except ValueError:
                        thresholds[hmm] = {'threshold': 50.0, 'type': 'full'}
    return thresholds


def load_hmm_template_thresholds(hmm_template_file: Path) -> dict[str, dict]:
    """Load custom HMM thresholds from hmm_table_template.txt."""
    thresholds = {}
    with open(hmm_template_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 11:
                hmm_col = parts[5] if len(parts) > 5 else ""
                threshold_col = parts[10] if len(parts) > 10 else ""
                
                if not hmm_col or re.match(r'^K\d{5}\.hmm$', hmm_col):
                    continue
                
                if '|' in threshold_col:
                    t_val, t_type = threshold_col.split('|', 1)
                    try:
                        thresholds[hmm_col] = {'threshold': float(t_val), 'type': t_type}
                    except ValueError:
                        thresholds[hmm_col] = {'threshold': 50.0, 'type': 'full'}
                elif threshold_col:
                    try:
                        thresholds[hmm_col] = {'threshold': float(threshold_col), 'type': 'full'}
                    except ValueError:
                        pass
    return thresholds


def load_motifs(motif_file: Optional[Path]) -> dict[str, re.Pattern]:
    """Load motif regex definitions."""
    motifs = {}
    if not motif_file:
        return motifs
    
    with open(motif_file) as f:
        for line in f:
            parts = line.strip().split(":")
            if len(parts) == 2:
                hmm, pattern = parts
                pattern = pattern.replace("X", "[ARNDCQEGHILKMFPSTWYV]")
                motifs[hmm] = re.compile(pattern)
    return motifs


def load_motif_pairs(pair_file: Optional[Path]) -> dict[str, str]:
    """Load motif pairs."""
    pairs = {}
    if not pair_file:
        return pairs
    
    with open(pair_file) as f:
        for line in f:
            parts = line.strip().split(":")
            if len(parts) == 2:
                pairs[parts[0]] = parts[1]
    return pairs


def parse_tblout(tblout_file: Path) -> pl.DataFrame:
    """Parse hmmsearch tblout format using Polars."""
    data = []
    with open(tblout_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 10:
                continue
            data.append({
                'seq_id': parts[0],
                'hmm_name': parts[2],
                'full_score': float(parts[5]),
                'domain_score': float(parts[8])
            })
    return pl.DataFrame(data) if data else pl.DataFrame(schema={
        'seq_id': pl.Utf8, 'hmm_name': pl.Utf8, 
        'full_score': pl.Float64, 'domain_score': pl.Float64
    })


@app.command()
def main(
    kofam_results: Path = typer.Option(..., "--kofam_results", help="KOfam hmmsearch tblout"),
    proteins: Path = typer.Option(..., "--proteins", help="Fasta file for motif validation"),
    ko_list: Path = typer.Option(..., "--ko_list", help="KOfam ko_list file with thresholds"),
    hmm_template: Path = typer.Option(..., "--hmm_template", help="hmm_table_template.txt"),
    contig_map: Path = typer.Option(..., "--contig_map", help="contigID to genomeID mapping TSV"),
    output: Path = typer.Option(..., "--output", help="Output TSV file"),
    custom_results: Optional[Path] = typer.Option(None, "--custom_results", help="Custom hmmsearch tblout"),
    motif_file: Optional[Path] = typer.Option(None, "--motif_file", help="Motif definitions"),
    motif_pair_file: Optional[Path] = typer.Option(None, "--motif_pair_file", help="Motif pairs"),
    cluster_tsv: Optional[Path] = typer.Option(None, "--cluster_tsv", help="MMseqs2 cluster TSV"),
):
    """Parse HMMsearch results for METABOLIC pipeline."""
    
    console.log("[bold blue]Loading reference data...")
    
    # 1. Load mappings
    contig_to_genome = load_contig_map(contig_map)
    thresholds = load_ko_list_thresholds(ko_list)
    custom_thresholds = load_hmm_template_thresholds(hmm_template)
    thresholds.update(custom_thresholds)
    
    motif_regex = load_motifs(motif_file)
    motif_pairs = load_motif_pairs(motif_pair_file)
    
    console.log(f"[green]Loaded {len(thresholds)} HMM thresholds")
    
    # 2. Parse Results
    console.log("[bold blue]Parsing HMM results...")
    df = parse_tblout(kofam_results)
    
    if custom_results:
        df_cust = parse_tblout(custom_results)
        df = pl.concat([df, df_cust])
    
    if df.height == 0:
        console.log("[yellow]No hits found, writing empty output")
        pl.DataFrame(schema={
            'seq_id': pl.Utf8, 'hmm_name': pl.Utf8, 
            'full_score': pl.Float64, 'domain_score': pl.Float64, 'genome_id': pl.Utf8
        }).write_csv(output, separator='\t')
        return

    console.log(f"[green]Parsed {df.height} raw hits")
    
    # 3. Pre-load sequences for motif validation
    seq_dict = SeqIO.to_dict(SeqIO.parse(proteins, "fasta"))
    
    # 4. Filter by Thresholds and Motifs (row-wise logic required)
    df_dict = df.to_dicts()
    score_map = {(r['seq_id'], r['hmm_name']): r['full_score'] for r in df_dict}
    
    valid_hits = []
    for row in df_dict:
        hmm = f"{row['hmm_name']}.hmm"
        
        if hmm not in thresholds:
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
            
        # Motif validation
        if row['hmm_name'] in motif_regex:
            seq_record = seq_dict.get(row['seq_id'])
            if seq_record:
                seq_str = str(seq_record.seq)
                if motif_regex[row['hmm_name']].search(seq_str):
                    valid_hits.append(row)
        elif row['hmm_name'] in motif_pairs:
            partner_hmm = motif_pairs[row['hmm_name']]
            my_score = row['full_score']
            partner_score = score_map.get((row['seq_id'], partner_hmm), 0.0)
            if my_score >= partner_score:
                valid_hits.append(row)
        else:
            valid_hits.append(row)

    console.log(f"[green]Filtered to {len(valid_hits)} valid hits")
    
    # 5. Expand hits if clustering was used
    if cluster_tsv:
        cluster_map = load_cluster_map(cluster_tsv)
        expanded_rows = []
        for row in valid_hits:
            rep_id = row['seq_id']
            members = cluster_map.get(rep_id, [rep_id])
            for member_id in members:
                new_row = row.copy()
                new_row['seq_id'] = member_id
                expanded_rows.append(new_row)
        valid_hits = expanded_rows
        console.log(f"[green]Expanded to {len(valid_hits)} hits after cluster expansion")
    
    # 6. Add genome_id and write output
    if valid_hits:
        final_df = pl.DataFrame(valid_hits)
        final_df = final_df.with_columns(
            pl.col('seq_id').map_elements(
                extract_contig_from_protein, return_dtype=pl.Utf8
            ).alias('contig_id')
        ).with_columns(
            pl.col('contig_id').replace(contig_to_genome, default='UNKNOWN').alias('genome_id')
        ).drop('contig_id')
        
        final_df.write_csv(output, separator='\t')
        console.log(f"[bold green]✓ Wrote {final_df.height} hits to {output}")
    else:
        pl.DataFrame(schema={
            'seq_id': pl.Utf8, 'hmm_name': pl.Utf8, 
            'full_score': pl.Float64, 'domain_score': pl.Float64, 'genome_id': pl.Utf8
        }).write_csv(output, separator='\t')
        console.log("[yellow]No valid hits, wrote empty output")


if __name__ == "__main__":
    app()
