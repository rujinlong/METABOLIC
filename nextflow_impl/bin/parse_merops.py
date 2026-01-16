#!/usr/bin/env python3
"""
Parse Diamond BLASTP results against MEROPS database.
Maps protein hits to peptidase families using contig_map for genome resolution.
Supports cluster expansion from MMseqs2 clustering.
"""

import re
from pathlib import Path
from typing import Optional

import polars as pl
import typer
from rich.console import Console

console = Console(stderr=True)
app = typer.Typer(add_completion=False)


def load_merops_mapping(lib_file: Path) -> dict[str, str]:
    """Load pepunit.lib to map Sequence ID -> Family."""
    mapping = {}
    try:
        with open(lib_file, 'r', errors='replace') as f:
            for line in f:
                if line.startswith(">"):
                    line = line.strip().replace('\r', '')
                    parts = line.split(maxsplit=1)
                    if not parts:
                        continue
                    full_id = parts[0][1:]  # remove >
                    
                    match = re.search(r'\#(.+?)\#', line)
                    if match:
                        family = match.group(1)
                        mapping[full_id] = family
    except FileNotFoundError:
        console.log(f"[red]Warning: {lib_file} not found")
    return mapping


def extract_contig_from_protein(protein_id: str) -> str:
    """Extract contigID from proteinID (format: contigID_proteinNumber)."""
    parts = protein_id.rsplit('_', 1)
    return parts[0] if len(parts) > 1 else protein_id


def load_contig_map(contig_map_file: Path) -> dict[str, str]:
    """Load contig-to-genome mapping."""
    contig_to_genome = {}
    with open(contig_map_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                contig_to_genome[parts[0]] = parts[1]
    return contig_to_genome


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


@app.command()
def main(
    input_file: Path = typer.Argument(..., help="Diamond m8 output file"),
    lib_file: Path = typer.Argument(..., help="pepunit.lib file"),
    output: Path = typer.Argument(..., help="Output TSV file"),
    contig_map: Path = typer.Option(..., "--contig_map", help="contigID to genomeID mapping TSV"),
    cluster_tsv: Optional[Path] = typer.Option(None, "--cluster_tsv", help="MMseqs2 cluster TSV"),
):
    """Parse Diamond BLASTP results against MEROPS database."""
    
    console.log("[bold blue]Loading MEROPS mapping...")
    merops_mapping = load_merops_mapping(lib_file)
    console.log(f"[green]Loaded {len(merops_mapping)} MEROPS families")
    
    contig_to_genome = load_contig_map(contig_map)
    cluster_map = load_cluster_map(cluster_tsv) if cluster_tsv else None
    
    # Parse m8 file
    console.log("[bold blue]Parsing Diamond results...")
    rows = []
    
    with open(input_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            
            seq_id = parts[0]
            target_id = parts[1]
            evalue = parts[10]
            
            family = merops_mapping.get(target_id)
            if not family:
                continue
            
            # Expand to all cluster members if clustering was used
            seq_ids = cluster_map.get(seq_id, [seq_id]) if cluster_map else [seq_id]
            
            for member_id in seq_ids:
                contig_id = extract_contig_from_protein(member_id)
                genome_id = contig_to_genome.get(contig_id, 'UNKNOWN')
                rows.append({
                    'seq_id': member_id,
                    'family': family,
                    'target_id': target_id,
                    'evalue': evalue,
                    'genome_id': genome_id
                })
    
    # Write output
    if rows:
        df = pl.DataFrame(rows)
        df.write_csv(output, separator='\t')
        console.log(f"[bold green]✓ Wrote {df.height} hits to {output}")
    else:
        pl.DataFrame(schema={
            'seq_id': pl.Utf8, 'family': pl.Utf8, 'target_id': pl.Utf8,
            'evalue': pl.Utf8, 'genome_id': pl.Utf8
        }).write_csv(output, separator='\t')
        console.log("[yellow]No hits found, wrote empty output")


if __name__ == "__main__":
    app()
