#!/usr/bin/env python3
"""
Parse HMMsearch results for METABOLIC.

This script filters HMM hits based on:
1. Per-HMM thresholds (from ko_list and hmm_table_template)
2. Motif validation (active site residues)
3. Motif pair comparisons (competing HMMs)
4. Custom DB priority (prefer custom HMM hits over KOfam for same function)

It uses a contig_map to resolve protein IDs to genome/MAG IDs.
"""
from __future__ import annotations

import sys
import re
from pathlib import Path
from typing import Optional

import polars as pl
import typer
from Bio import SeqIO
from rich.console import Console

console = Console(stderr=True)
app = typer.Typer(add_completion=False)


def load_contig_map(contig_map_file: Path) -> dict[str, str]:
    """Load contig-to-genome mapping using Polars."""
    df = pl.read_csv(
        contig_map_file, 
        separator='\t', 
        has_header=False, 
        new_columns=['contig_id', 'genome_id']
    )
    return dict(zip(df['contig_id'].to_list(), df['genome_id'].to_list()))


def extract_contig_from_protein(protein_id: str) -> str:
    """Extract contigID from proteinID (format: contigID_proteinNumber)."""
    parts = protein_id.rsplit('_', 1)
    return parts[0] if len(parts) > 1 else protein_id


def load_cluster_map(cluster_tsv_file: Path) -> dict[str, list[str]]:
    """Load MMseqs2 cluster membership using Polars: rep_id -> [member_ids]."""
    df = pl.read_csv(
        cluster_tsv_file, 
        separator='\t', 
        has_header=False, 
        new_columns=['rep_id', 'member_id']
    )
    # Group by rep_id and aggregate member_ids into list
    grouped = df.group_by('rep_id').agg(pl.col('member_id').alias('members'))
    return {row['rep_id']: row['members'] for row in grouped.iter_rows(named=True)}


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


def load_hmm_template_thresholds(hmm_template_file: Path) -> tuple[dict[str, dict], dict[str, str]]:
    """
    Load custom HMM thresholds and HMM-to-KO mapping from hmm_table_template.txt.
    Returns: (thresholds_dict, hmm_to_ko_dict)
    
    Column 6 (index 5) = HMM filename
    Column 7 (index 6) = Corresponding KO
    Column 11 (index 10) = threshold|score_type
    """
    thresholds = {}
    hmm_to_ko = {}  # custom HMM -> corresponding KO
    
    with open(hmm_template_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 11:
                hmm_col = parts[5] if len(parts) > 5 else ""
                ko_col = parts[6] if len(parts) > 6 else ""
                threshold_col = parts[10] if len(parts) > 10 else ""
                
                # Skip empty or KOfam entries (we only want custom HMMs here)
                if not hmm_col or re.match(r'^K\d{5}\.hmm$', hmm_col):
                    continue
                
                # Build HMM -> KO mapping for deduplication
                if ko_col and re.match(r'^K\d{5}$', ko_col):
                    hmm_to_ko[hmm_col] = ko_col
                
                # Parse threshold
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
    
    return thresholds, hmm_to_ko


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


def parse_tblout(tblout_file: Path, source: str) -> pl.DataFrame:
    """Parse hmmsearch tblout format using Polars, with source tracking."""
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
                'domain_score': float(parts[8]),
                'source': source  # 'kofam' or 'custom'
            })
    return pl.DataFrame(data) if data else pl.DataFrame(schema={
        'seq_id': pl.Utf8, 'hmm_name': pl.Utf8, 
        'full_score': pl.Float64, 'domain_score': pl.Float64, 'source': pl.Utf8
    })


def deduplicate_hits(
    valid_hits: list[dict], 
    hmm_to_ko: dict[str, str],
    prefer_custom: bool
) -> list[dict]:
    """
    Deduplicate hits when same protein matches both custom and KOfam for same function.
    
    If prefer_custom=True:
    - When a protein has hits from both custom_db and KOfam for the same underlying KO,
      keep only the custom_db hit.
    """
    if not prefer_custom or not hmm_to_ko:
        return valid_hits
    
    # Build reverse map: KO -> [custom HMMs that map to it]
    ko_to_custom_hmms = {}
    for hmm, ko in hmm_to_ko.items():
        if ko not in ko_to_custom_hmms:
            ko_to_custom_hmms[ko] = set()
        ko_to_custom_hmms[ko].add(hmm.replace('.hmm', ''))
    
    # Group hits by seq_id
    hits_by_seq: dict[str, list[dict]] = {}
    for hit in valid_hits:
        seq_id = hit['seq_id']
        if seq_id not in hits_by_seq:
            hits_by_seq[seq_id] = []
        hits_by_seq[seq_id].append(hit)
    
    # For each seq_id, check for overlapping KOfam/custom hits
    deduplicated = []
    removed_count = 0
    
    for seq_id, hits in hits_by_seq.items():
        # Identify KOfam hits and custom hits
        kofam_hits = [h for h in hits if h.get('source') == 'kofam']
        custom_hits = [h for h in hits if h.get('source') == 'custom']
        
        # For each custom hit, find if there's a corresponding KOfam hit
        custom_kos = set()
        for ch in custom_hits:
            hmm_name = ch['hmm_name']
            hmm_file = f"{hmm_name}.hmm"
            if hmm_file in hmm_to_ko:
                custom_kos.add(hmm_to_ko[hmm_file])
        
        # Filter KOfam hits: remove if custom already covers that KO
        filtered_kofam = []
        for kh in kofam_hits:
            ko = kh['hmm_name']  # KOfam hmm_name is the KO id (e.g., K00001)
            if ko in custom_kos:
                removed_count += 1
                continue  # Skip this KOfam hit, custom takes priority
            filtered_kofam.append(kh)
        
        deduplicated.extend(custom_hits)
        deduplicated.extend(filtered_kofam)
    
    if removed_count > 0:
        console.log(f"[cyan]Custom DB priority: removed {removed_count} redundant KOfam hits")
    
    return deduplicated


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
    prefer_custom: bool = typer.Option(True, "--prefer_custom/--no-prefer-custom", 
                                       help="Prefer custom_db hits over KOfam when both match same function"),
):
    """Parse HMMsearch results for METABOLIC pipeline."""
    
    console.log("[bold blue]Loading reference data...")
    
    # 1. Load mappings
    contig_to_genome = load_contig_map(contig_map)
    thresholds = load_ko_list_thresholds(ko_list)
    custom_thresholds, hmm_to_ko = load_hmm_template_thresholds(hmm_template)
    thresholds.update(custom_thresholds)
    
    motif_regex = load_motifs(motif_file)
    motif_pairs = load_motif_pairs(motif_pair_file)
    
    console.log(f"[green]Loaded {len(thresholds)} HMM thresholds, {len(hmm_to_ko)} custom-to-KO mappings")
    
    # 2. Parse Results with source tracking
    console.log("[bold blue]Parsing HMM results...")
    df = parse_tblout(kofam_results, source='kofam')
    
    if custom_results:
        df_cust = parse_tblout(custom_results, source='custom')
        df = pl.concat([df, df_cust])
    
    if df.height == 0:
        console.log("[yellow]No hits found, writing empty output")
        pl.DataFrame(schema={
            'seq_id': pl.Utf8, 'hmm_name': pl.Utf8, 
            'full_score': pl.Float64, 'domain_score': pl.Float64, 
            'source': pl.Utf8, 'genome_id': pl.Utf8
        }).write_csv(output, separator='\t')
        return

    console.log(f"[green]Parsed {df.height} raw hits (KOfam + Custom)")
    
    # 3. Pre-load sequences for motif validation
    seq_dict = SeqIO.to_dict(SeqIO.parse(proteins, "fasta"))
    
    # 4. Filter by Thresholds and Motifs
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
    
    # 5. Apply custom DB priority deduplication
    if prefer_custom and custom_results:
        valid_hits = deduplicate_hits(valid_hits, hmm_to_ko, prefer_custom)
        console.log(f"[green]After deduplication: {len(valid_hits)} hits")
    
    # 6. Expand hits if clustering was used
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
    
    # 7. Add genome_id and write output
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
            'full_score': pl.Float64, 'domain_score': pl.Float64,
            'source': pl.Utf8, 'genome_id': pl.Utf8
        }).write_csv(output, separator='\t')
        console.log("[yellow]No valid hits, wrote empty output")


if __name__ == "__main__":
    app()
