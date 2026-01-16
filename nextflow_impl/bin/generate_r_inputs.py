#!/usr/bin/env python3
"""
Generate R_input files for METABOLIC biogeochemical cycle visualization.

This script creates per-genome R_input.txt files that are used by
draw_biogeochemical_cycles.R to generate N/C/S/Other cycle PDFs.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import polars as pl
import typer
from rich.console import Console

console = Console(stderr=True)
app = typer.Typer(add_completion=False)


def load_r_pathways(r_pathways_file: Path) -> dict[str, dict]:
    """
    Load R_pathways.txt defining biogeochemical cycle steps.
    
    Format: STEP_ID:STEP_NAME<TAB>HMM1,HMM2;HMM3  (semicolon = AND, NO| = NOT)
    
    Returns: dict[step_id] = {'name': str, 'hmm_groups': list of (hmms, is_not)}
    """
    pathways = {}
    with open(r_pathways_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            
            step_info = parts[0]
            hmm_string = parts[1]
            
            # Parse step ID and name
            step_id, step_name = step_info.split(':', 1) if ':' in step_info else (step_info, step_info)
            
            # Parse HMM groups (separated by ;)
            hmm_groups = []
            for group in hmm_string.split(';'):
                group = group.strip()
                if not group:
                    continue
                # Check for NOT logic (NO|)
                hmms_in_group = []
                for hmm in group.split(','):
                    hmm = hmm.strip()
                    if hmm.startswith('NO|'):
                        hmms_in_group.append((hmm[3:], True))  # (hmm_name, is_not)
                    else:
                        hmms_in_group.append((hmm, False))
                hmm_groups.append(hmms_in_group)
            
            pathways[step_id] = {
                'name': step_name,
                'hmm_groups': hmm_groups
            }
    
    return pathways


def check_pathway_presence(
    hmm_groups: list,
    genome_hmms: set[str]
) -> bool:
    """
    Check if a genome has the required HMMs for a pathway.
    
    Logic:
    - Groups are ANDed together (all groups must pass)
    - Within a group, HMMs are ORed (any HMM match is enough)
    - NO|HMM means the HMM must NOT be present
    """
    if not hmm_groups:
        return False
    
    for group in hmm_groups:
        # Within group: check if any positive HMM matches,
        # and no negative HMM matches
        group_pass = False
        has_positive = False
        
        for hmm, is_not in group:
            hmm_base = hmm.replace('.hmm', '')
            hmm_with_ext = f"{hmm_base}.hmm" if not hmm.endswith('.hmm') else hmm
            
            # Check both with and without .hmm extension
            hmm_present = hmm_base in genome_hmms or hmm_with_ext in genome_hmms
            
            if is_not:
                # NOT logic: if present, group fails
                if hmm_present:
                    group_pass = False
                    break
            else:
                has_positive = True
                if hmm_present:
                    group_pass = True
        
        # If no positive match in this AND-group, pathway fails
        if has_positive and not group_pass:
            return False
    
    return True


def load_genome_hmms(hits_tsv: Path) -> dict[str, set[str]]:
    """Load HMM hits per genome."""
    genome_hmms: dict[str, set[str]] = {}
    
    df = pl.read_csv(hits_tsv, separator='\t')
    
    for row in df.iter_rows(named=True):
        genome_id = row.get('genome_id', 'UNKNOWN')
        hmm_name = row.get('hmm_name', '')
        
        if genome_id not in genome_hmms:
            genome_hmms[genome_id] = set()
        
        # Store both with and without .hmm extension
        genome_hmms[genome_id].add(hmm_name)
        if not hmm_name.endswith('.hmm'):
            genome_hmms[genome_id].add(f"{hmm_name}.hmm")
    
    return genome_hmms


@app.command()
def main(
    hits_tsv: Path = typer.Argument(..., help="Parsed HMM hits TSV (all_hits_parsed.tsv)"),
    r_pathways: Path = typer.Argument(..., help="R_pathways.txt file"),
    outdir: Path = typer.Option("Nutrient_Cycling_Diagram_Input", "--outdir", help="Output directory"),
):
    """Generate R_input files for biogeochemical cycle visualization."""
    
    console.log("[bold blue]Loading data...")
    
    # Load pathways and HMM hits
    pathways = load_r_pathways(r_pathways)
    genome_hmms = load_genome_hmms(hits_tsv)
    
    console.log(f"[green]Loaded {len(pathways)} pathways, {len(genome_hmms)} genomes")
    
    # Create output directory
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Track totals for summary file
    pathway_counts: dict[str, int] = {step_id: 0 for step_id in pathways}
    total_genomes = len(genome_hmms)
    
    # Generate per-genome R_input files
    for genome_id, hmms in genome_hmms.items():
        if genome_id == 'UNKNOWN':
            continue
        
        output_file = outdir / f"{genome_id}.R_input.txt"
        
        with open(output_file, 'w') as f:
            for step_id in sorted(pathways.keys()):
                pathway = pathways[step_id]
                is_present = check_pathway_presence(pathway['hmm_groups'], hmms)
                presence_val = 1 if is_present else 0
                
                f.write(f"{step_id}:{pathway['name']}\t{presence_val}\n")
                
                if is_present:
                    pathway_counts[step_id] += 1
        
        console.log(f"[dim]Generated: {output_file.name}")
    
    # Generate Total.R_input.txt
    total_file = outdir / "Total.R_input.txt"
    with open(total_file, 'w') as f:
        for step_id in sorted(pathways.keys()):
            pathway = pathways[step_id]
            count = pathway_counts[step_id]
            percentage = count / total_genomes if total_genomes > 0 else 0
            
            f.write(f"{step_id}:{pathway['name']}\t{count}\t{percentage:.4f}\n")
    
    console.log(f"[bold green]✓ Generated {len(genome_hmms)} R_input files + Total.R_input.txt")


if __name__ == "__main__":
    app()
