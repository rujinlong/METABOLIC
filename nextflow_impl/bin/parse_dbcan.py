#!/usr/bin/env python3
"""
Parse dbCAN hmmscan domtblout results for METABOLIC.
Implements the legacy hmmscan-parser-dbCANmeta.py logic with overlap resolution.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import polars as pl
import typer
from rich.console import Console

console = Console(stderr=True)
app = typer.Typer(add_completion=False)

# Legacy defaults from hmmscan-parser-dbCANmeta.py
DEFAULT_EVAL = 1e-15
DEFAULT_COVERAGE = 0.35


def extract_contig_from_protein(protein_id: str) -> str:
    """Extract contigID from proteinID (format: contigID_proteinNumber)."""
    parts = protein_id.rsplit('_', 1)
    return parts[0] if len(parts) > 1 else protein_id


def load_contig_map(contig_map_file: Path) -> dict[str, str]:
    """Load contig-to-genome mapping using Polars."""
    df = pl.read_csv(
        contig_map_file, 
        separator='\t', 
        has_header=False, 
        new_columns=['contig_id', 'genome_id']
    )
    # plot head of df
    console.log(df.head())
    return dict(zip(df['contig_id'].to_list(), df['genome_id'].to_list()))


def load_cluster_map(cluster_tsv_file: Path) -> dict[str, list[str]]:
    """Load MMseqs2 cluster membership using Polars: rep_id -> [member_ids]."""
    df = pl.read_csv(
        cluster_tsv_file, 
        separator='\t', 
        has_header=False, 
        new_columns=['rep_id', 'member_id']
    )
    grouped = df.group_by('rep_id').agg(pl.col('member_id').alias('members'))
    return {row['rep_id']: row['members'] for row in grouped.iter_rows(named=True)}


def parse_domtblout(input_file: Path) -> list[dict]:
    """
    Parse hmmscan domtblout format.
    
    In hmmscan output:
    - Column 0 (target): HMM name (e.g., AA3.hmm)
    - Column 3 (query): Protein sequence ID (e.g., BBKHAN_03036)
    """
    hits = []
    with open(input_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            
            hits.append({
                'hmm_name': parts[0],      # HMM name (target in hmmscan)
                'hmm_len': int(parts[2]),  # HMM length
                'seq_id': parts[3],        # Protein ID (query in hmmscan)
                'seq_len': int(parts[5]),  # Protein length
                'evalue': float(parts[12]),  # i-Evalue
                'hmm_from': int(parts[15]),
                'hmm_to': int(parts[16]),
                'ali_from': int(parts[17]),
                'ali_to': int(parts[18])
            })
    return hits


def resolve_overlaps(hits_by_hmm: dict[str, list[dict]]) -> list[dict]:
    """Resolve overlapping domain hits within each HMM group."""
    final_hits = []
    
    for hmm, hits in hits_by_hmm.items():
        # Sort by alignment coordinates
        hits.sort(key=lambda x: (x['ali_from'], x['ali_to']))
        
        i = 0
        while i < len(hits) - 1:
            curr = hits[i]
            next_h = hits[i + 1]
            
            len1 = curr['ali_to'] - curr['ali_from']
            len2 = next_h['ali_to'] - next_h['ali_from']
            
            overlap_end = min(curr['ali_to'], next_h['ali_to'])
            overlap_start = next_h['ali_from']
            len3 = overlap_end - overlap_start
            
            if len3 > 0:
                ratio1 = len3 / len1 if len1 > 0 else 0
                ratio2 = len3 / len2 if len2 > 0 else 0
                
                if ratio1 > 0.5 or ratio2 > 0.5:
                    # Remove hit with worse E-value
                    if curr['evalue'] < next_h['evalue']:
                        hits.pop(i + 1)
                    else:
                        hits.pop(i)
                        if i > 0:
                            i -= 1
                    continue
            i += 1
            
        final_hits.extend(hits)
    
    return final_hits


@app.command()
def main(
    input_file: Path = typer.Argument(..., help="hmmscan domtblout file"),
    output: Path = typer.Argument(..., help="Output TSV file"),
    contig_map: Path = typer.Option(..., "--contig_map", help="contigID to genomeID mapping TSV"),
    cluster_tsv: Optional[Path] = typer.Option(None, "--cluster_tsv", help="MMseqs2 cluster TSV"),
    evalue_cutoff: float = typer.Option(DEFAULT_EVAL, "--evalue", help="E-value cutoff"),
    coverage_cutoff: float = typer.Option(DEFAULT_COVERAGE, "--coverage", help="Coverage cutoff"),
):
    """Parse dbCAN hmmscan results for METABOLIC pipeline."""
    
    console.log("[bold blue]Parsing dbCAN domtblout...")
    
    # Load mappings
    contig_to_genome = load_contig_map(contig_map)
    cluster_map = load_cluster_map(cluster_tsv) if cluster_tsv else None
    
    # Parse and group by HMM
    raw_hits = parse_domtblout(input_file)
    console.log(f"[green]Parsed {len(raw_hits)} raw domain hits")
    
    # Group by HMM name for overlap resolution
    hits_by_hmm: dict[str, list[dict]] = {}
    for hit in raw_hits:
        if hit['hmm_name'] not in hits_by_hmm:
            hits_by_hmm[hit['hmm_name']] = []
        hits_by_hmm[hit['hmm_name']].append(hit)
    
    # Resolve overlaps
    final_hits = resolve_overlaps(hits_by_hmm)
    console.log(f"[green]After overlap resolution: {len(final_hits)} hits")
    
    # Filter by E-value and coverage, expand clusters, collect rows
    rows = []
    for h in final_hits:
        cov = (h['hmm_to'] - h['hmm_from']) / h['hmm_len'] if h['hmm_len'] > 0 else 0
        
        if h['evalue'] <= evalue_cutoff and cov >= coverage_cutoff:
            # Use seq_id (protein ID) for genome resolution, NOT hmm_name!
            rep_id = h['seq_id']
            seq_ids = cluster_map.get(rep_id, [rep_id]) if cluster_map else [rep_id]
            
            for seq_id in seq_ids:
                contig_id = extract_contig_from_protein(seq_id)
                genome_id = contig_to_genome.get(contig_id, 'UNKNOWN')
                rows.append({
                    'seq_id': seq_id,
                    'hmm_name': h['hmm_name'].replace('.hmm', ''),  # Strip .hmm suffix
                    'hmm_len': h['hmm_len'],
                    'seq_len': h['seq_len'],
                    'evalue': h['evalue'],
                    'hmm_from': h['hmm_from'],
                    'hmm_to': h['hmm_to'],
                    'ali_from': h['ali_from'],
                    'ali_to': h['ali_to'],
                    'coverage': round(cov, 4),
                    'genome_id': genome_id
                })
    
    # Write output
    if rows:
        df = pl.DataFrame(rows)
        df.write_csv(output, separator='\t')
        console.log(f"[bold green]✓ Wrote {df.height} hits to {output}")
    else:
        pl.DataFrame(schema={
            'seq_id': pl.Utf8, 'hmm_name': pl.Utf8, 'hmm_len': pl.Int64, 'seq_len': pl.Int64,
            'evalue': pl.Float64, 'hmm_from': pl.Int64, 'hmm_to': pl.Int64,
            'ali_from': pl.Int64, 'ali_to': pl.Int64, 'coverage': pl.Float64, 'genome_id': pl.Utf8
        }).write_csv(output, separator='\t')
        console.log("[yellow]No hits passed filters, wrote empty output")


if __name__ == "__main__":
    app()
