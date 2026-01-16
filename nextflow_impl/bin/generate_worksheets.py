#!/usr/bin/env python3
"""
Generate METABOLIC result worksheets.
Produces 6 worksheets summarizing HMM hits, KEGG modules, dbCAN, and MEROPS results.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

import polars as pl
import typer
from rich.console import Console
from rich.progress import track

console = Console(stderr=True)
app = typer.Typer(add_completion=False)


def load_hmm_template(template_file: Path) -> dict:
    """Load hmm_table_template.txt."""
    templates = {}
    with open(template_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            entry = parts[0]
            templates[entry] = {
                'category': parts[1],
                'function': parts[2],
                'gene_abbr': parts[3],
                'gene_name': parts[4],
                'hmm_file': parts[5],
                'ko': parts[6],
                'reaction': parts[7],
                'product': parts[9],
                'threshold': parts[10] if len(parts) > 10 else ""
            }
    return templates


def load_hmm_template_2(template_file: Path) -> list[dict]:
    """Load hmm_table_template_2.txt."""
    templates = []
    with open(template_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            templates.append({
                'id': parts[0],
                'entry_refs': parts[1],
                'category': parts[2],
                'function': parts[3],
                'gene_abbr': parts[4]
            })
    return templates


def load_kegg_module_db(module_file: Path) -> dict[str, str]:
    """Load ko00002.keg to get Module -> Category mapping."""
    module_cat = {}
    current_cat = ""
    with open(module_file) as f:
        for line in f:
            if line.startswith("C"):
                current_cat = line.strip()[1:].strip()
            elif line.startswith("D"):
                parts = line.strip()[1:].strip().split()
                if parts:
                    module_cat[parts[0]] = current_cat
    return module_cat


def load_kegg_steps(step_file: Path) -> dict:
    """Load kegg_module_step_db.txt."""
    steps = {}
    with open(step_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            desc = parts[0]
            k_string = parts[1]
            step_id = parts[2]
            
            # Parse step_id: M00001+01 -> module=M00001
            if "+" in step_id:
                module_id = step_id.split("+")[0]
            else:
                module_id = step_id
                
            steps[step_id] = {
                'desc': desc,
                'k_string': k_string,
                'module': module_id
            }
    return steps


def check_step_presence(k_string: str, present_kos: set[str]) -> bool:
    """Evaluate boolean KO string like '(K001 or K002) and K003'."""
    def replacer(match: re.Match) -> str:
        ko = match.group(0)
        return "True" if ko in present_kos else "False"
    
    expr = re.sub(r'K\d{5}', replacer, k_string)
    expr = expr.replace(" and ", " and ").replace(" or ", " or ")
    
    try:
        return eval(expr)
    except Exception:
        return False


@app.command()
def main(
    hits: Path = typer.Option(..., "--hits", help="Parsed HMM hits TSV"),
    template1: Path = typer.Option(..., "--template1", help="hmm_table_template.txt"),
    template2: Path = typer.Option(..., "--template2", help="hmm_table_template_2.txt"),
    kegg_module_db: Path = typer.Option(..., "--kegg_module_db", help="ko00002.keg"),
    kegg_steps: Path = typer.Option(..., "--kegg_steps", help="kegg_module_step_db.txt"),
    output_dir: Path = typer.Option(..., "--output_dir", help="Output directory"),
    genome_ids: str = typer.Option(..., "--genome_ids", help="Comma-separated genome IDs"),
    dbcan_hits: Optional[Path] = typer.Option(None, "--dbcan_hits", help="Parsed dbCAN hits TSV"),
    merops_hits: Optional[Path] = typer.Option(None, "--merops_hits", help="Parsed MEROPS hits TSV"),
):
    """Generate METABOLIC result worksheets."""
    
    genome_list = [g.strip() for g in genome_ids.split(",")]
    console.log(f"[bold blue]Generating worksheets for {len(genome_list)} genomes...")
    
    # Load HMM hits
    try:
        df_hits = pl.read_csv(hits, separator='\t')
    except Exception:
        df_hits = pl.DataFrame(schema={'seq_id': pl.Utf8, 'hmm_name': pl.Utf8, 'genome_id': pl.Utf8})
    
    console.log(f"[green]Loaded {df_hits.height} HMM hits")
    
    # Build presence and KO maps
    presence_map: dict[str, dict[str, int]] = {g: {} for g in genome_list}
    ko_map: dict[str, set[str]] = {g: set() for g in genome_list}
    seq_map: dict[str, str] = {}
    
    if df_hits.height > 0:
        for row in df_hits.iter_rows(named=True):
            g = row.get('genome_id', '')
            h = row.get('hmm_name', '')
            seq_id = row.get('seq_id', '')
            
            seq_map[seq_id] = g
            
            h_clean = f"{h}.hmm" if not h.endswith(".hmm") else h
            
            if g in presence_map:
                presence_map[g][h_clean] = presence_map[g].get(h_clean, 0) + 1
                
                if h.startswith("K") and len(h) >= 6 and h[1:6].isdigit():
                    ko_map[g].add(h.split(".")[0])
    
    # Load dbCAN hits
    dbcan_counts: dict[str, dict[str, int]] = {g: {} for g in genome_list}
    if dbcan_hits and dbcan_hits.exists():
        try:
            df_dbcan = pl.read_csv(dbcan_hits, separator='\t')
            for row in df_dbcan.iter_rows(named=True):
                gid = row.get('genome_id', 'UNKNOWN')
                family = row.get('hmm_name', '')  # CAZyme family (e.g., GH1, AA3)
                if gid in dbcan_counts and family:
                    dbcan_counts[gid][family] = dbcan_counts[gid].get(family, 0) + 1
        except Exception as e:
            console.log(f"[yellow]Warning: Failed to parse dbCAN hits: {e}")
    
    # Load MEROPS hits
    merops_counts: dict[str, dict[str, int]] = {g: {} for g in genome_list}
    if merops_hits and merops_hits.exists():
        try:
            df_merops = pl.read_csv(merops_hits, separator='\t')
            for row in df_merops.iter_rows(named=True):
                gid = row.get('genome_id', 'UNKNOWN')
                family = row.get('family', '')
                if gid in merops_counts:
                    merops_counts[gid][family] = merops_counts[gid].get(family, 0) + 1
        except Exception as e:
            console.log(f"[yellow]Warning: Failed to parse MEROPS hits: {e}")
    
    # --- Worksheet 1 & 2 ---
    console.log("[bold blue]Generating Worksheet 1 & 2...")
    templates1 = load_hmm_template(template1)
    
    rows1 = []
    for entry_id, t in templates1.items():
        row = {
            'category': t['category'],
            'function': t['function'],
            'gene_abbr': t['gene_abbr'],
            'gene_name': t['gene_name'],
            'hmm_file': t['hmm_file'],
            'corresponding_ko': t['ko'],
            'reaction': t['reaction'],
            'product': t['product'],
            'hmm_threshold': t['threshold']
        }
        hmms = [x.strip() for x in t['hmm_file'].split(",")]
        for g in genome_list:
            count = sum(presence_map[g].get(h, 0) for h in hmms)
            row[f"{g}_presence"] = "Present" if count > 0 else "Absent"
            row[f"{g}_hit_count"] = count
        rows1.append(row)
    
    pl.DataFrame(rows1).write_csv(output_dir / "METABOLIC_result_worksheet1.tsv", separator='\t')
    
    templates2 = load_hmm_template_2(template2)
    rows2 = []
    for t in templates2:
        row = {
            'category': t['category'],
            'function': t['function'],
            'gene_abbr': t['gene_abbr']
        }
        refs = t['entry_refs'].split("||")
        target_hmms = []
        for r in refs:
            if r in templates1:
                t1_hmms = [x.strip() for x in templates1[r]['hmm_file'].split(",")]
                target_hmms.extend(t1_hmms)
        
        for g in genome_list:
            count = sum(presence_map[g].get(h, 0) for h in target_hmms)
            row[f"{g}_function_presence"] = "Present" if count > 0 else "Absent"
        rows2.append(row)
    
    pl.DataFrame(rows2).write_csv(output_dir / "METABOLIC_result_worksheet2.tsv", separator='\t')
    
    # --- Worksheet 4 (Step Presence) ---
    console.log("[bold blue]Generating Worksheet 4...")
    steps = load_kegg_steps(kegg_steps)
    module_cat = load_kegg_module_db(kegg_module_db)
    
    module_steps_count: dict[str, int] = defaultdict(int)
    module_steps_present: dict[str, dict[str, int]] = {m: {g: 0 for g in genome_list} for m in module_cat}
    
    rows4 = []
    for step_id in sorted(steps.keys()):
        step = steps[step_id]
        mod_id = step['module']
        cat = module_cat.get(mod_id, "")
        
        row = {
            'module_step': step_id,
            'module_name': step['desc'],
            'ko_expression': step['k_string'],
            'module_category': cat
        }
        
        module_steps_count[mod_id] += 1
        
        for g in genome_list:
            is_present = check_step_presence(step['k_string'], ko_map[g])
            row[f"{g}_step_presence"] = "Present" if is_present else "Absent"
            
            if is_present:
                if mod_id not in module_steps_present:
                    module_steps_present[mod_id] = {gid: 0 for gid in genome_list}
                module_steps_present[mod_id][g] += 1
        
        rows4.append(row)
    
    pl.DataFrame(rows4).write_csv(output_dir / "METABOLIC_result_worksheet4.tsv", separator='\t')
    
    # --- Worksheet 3 (Module Presence) ---
    console.log("[bold blue]Generating Worksheet 3...")
    all_modules = sorted(set(s['module'] for s in steps.values()))
    
    rows3 = []
    for mod_id in all_modules:
        name = ""
        for s_id in steps:
            if steps[s_id]['module'] == mod_id:
                name = steps[s_id]['desc']
                break
        
        cat = module_cat.get(mod_id, "")
        total_steps = module_steps_count[mod_id]
        
        row = {
            'module_id': mod_id,
            'module_name': name,
            'module_category': cat
        }
        
        for g in genome_list:
            present_count = module_steps_present.get(mod_id, {}).get(g, 0)
            if total_steps > 0 and present_count / total_steps >= 0.75:
                row[f"{g}_module_presence"] = "Present"
            else:
                row[f"{g}_module_presence"] = "Absent"
        rows3.append(row)
    
    pl.DataFrame(rows3).write_csv(output_dir / "METABOLIC_result_worksheet3.tsv", separator='\t')
    
    # --- Worksheet 5 (dbCAN / CAZyme) ---
    console.log("[bold blue]Generating Worksheet 5...")
    all_cazymes = set()
    for g in dbcan_counts:
        all_cazymes.update(dbcan_counts[g].keys())
    
    rows5 = []
    for fam in sorted(all_cazymes):
        row = {'cazyme_id': fam}
        for g in genome_list:
            row[f"{g}_hit_count"] = dbcan_counts[g].get(fam, 0)
        rows5.append(row)
    
    pl.DataFrame(rows5).write_csv(output_dir / "METABOLIC_result_worksheet5.tsv", separator='\t')
    
    # --- Worksheet 6 (MEROPS) ---
    console.log("[bold blue]Generating Worksheet 6...")
    all_pep = set()
    for g in merops_counts:
        all_pep.update(merops_counts[g].keys())
    
    rows6 = []
    for pep in sorted(all_pep):
        row = {'merops_id': pep}
        for g in genome_list:
            row[f"{g}_hit_count"] = merops_counts[g].get(pep, 0)
        rows6.append(row)
    
    pl.DataFrame(rows6).write_csv(output_dir / "METABOLIC_result_worksheet6.tsv", separator='\t')
    
    console.log("[bold green]✓ All worksheets generated successfully!")


if __name__ == "__main__":
    app()
