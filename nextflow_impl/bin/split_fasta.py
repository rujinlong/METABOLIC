#!/usr/bin/env python3
"""
Split a FASTA file into chunks of specified size.

Used for parallelizing HMM searches on large protein datasets.
"""
from __future__ import annotations

from pathlib import Path

import typer
from Bio import SeqIO
from rich.console import Console

console = Console(stderr=True)
app = typer.Typer(add_completion=False)


@app.command()
def main(
    input_fasta: Path = typer.Argument(..., help="Input FASTA file"),
    chunk_size: int = typer.Option(100000, "--chunk_size", help="Sequences per chunk"),
    prefix: str = typer.Option("chunk_", "--prefix", help="Output file prefix"),
):
    """Split FASTA into chunks for parallel processing."""
    
    console.log(f"[bold blue]Reading {input_fasta}...")
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    total_seqs = len(records)
    
    if total_seqs == 0:
        console.log("[red]No sequences found in input file!")
        raise typer.Exit(1)
    
    num_chunks = (total_seqs + chunk_size - 1) // chunk_size
    console.log(f"[green]Total sequences: {total_seqs}, splitting into {num_chunks} chunks of ≤{chunk_size}")
    
    for i in range(0, total_seqs, chunk_size):
        chunk_idx = i // chunk_size
        chunk = records[i:i + chunk_size]
        output_file = f"{prefix}{chunk_idx:04d}.faa"
        
        SeqIO.write(chunk, output_file, "fasta")
        console.log(f"[dim]  Wrote {len(chunk)} sequences to {output_file}")
    
    console.log(f"[bold green]✓ Created {num_chunks} chunk files")


if __name__ == "__main__":
    app()
