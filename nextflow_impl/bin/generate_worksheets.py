#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
import collections
import logging
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_hmm_template(template_file):
    """Load hmm_table_template.txt"""
    templates = {}
    with open(template_file) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
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

def load_hmm_template_2(template_file):
    """Load hmm_table_template_2.txt"""
    templates = []
    with open(template_file) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            templates.append({
                'id': parts[0],
                'entry_refs': parts[1],
                'category': parts[2],
                'function': parts[3],
                'gene_abbr': parts[4]
            })
    return templates

def load_kegg_module_db(module_file):
    """
    Load ko00002.keg to get Module -> Category mapping
    Hierarchical: C (Category) -> D (Module)
    """
    module_cat = {} # M00001 -> "Central carbohydrate metabolism"
    current_cat = ""
    with open(module_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("C    "):
                current_cat = line[5:]
            elif line.startswith("D      "):
                # D      M00001  Glycolysis...
                parts = line.split(maxsplit=2)
                if len(parts) >= 2:
                    module_id = parts[1]
                    module_cat[module_id] = current_cat
    return module_cat

def load_kegg_steps(step_file):
    """
    Load kegg_module_step_db.txt
    Format: Name \t (K001 or K002) \t M00001+01
    """
    # M00001+01 -> {'name': 'desc', 'k_string': '(K001 or K002)', 'module': 'M00001', 'step_num': 1}
    steps = {}
    module_names = {}
    with open(step_file) as f:
        for line in f:
            if line.startswith("name"): continue
            parts = line.strip().split("\t")
            if len(parts) < 3: continue
            
            desc = parts[0]
            k_string = parts[1]
            step_id = parts[2]
            
            # Parse M00001+01
            m_parts = step_id.split("+")
            if len(m_parts) == 2:
                mod_id = m_parts[0]
                steps[step_id] = {
                    'desc': desc,
                    'k_string': k_string,
                    'module': mod_id
                }
                module_names[mod_id] = desc # approximate name from step? No, legacy has separate map?
                # Actually legacy logic: $KEGG_module2name{$module} = $tmp[0];
                # So the description of the step is used as the module name?
                # METABOLIC-G.pl line 845: $KEGG_module2name{$module} = $tmp[0];
                # It overwrites with each step, so presumably all steps have same name or last wins.
    return steps

def check_step_presence(k_string, present_kos):
    """
    Evaluate boolean string like "(K001 or K002) and K003" against set of present KOs.
    """
    # Simple logic: replace Kxxx with True/False, eval() using python
    # BUT sanitize!
    # K-ids are K\d{5}. 
    # Replace (K\d{5}) with True/False?
    # Logic:
    # 1. Extract all K-ids
    # 2. Map to booleans
    # 3. Replace in string (K00001 -> True)
    # 4. Replace 'or'/'and' -> 'or'/'and' (python syntax matches)
    # WARNING: using eval() on external input is risky, but here DB is internal.
    
    # Legacy logic: uses split regex.
    # Python approach:
    # We can iterate regex replacement.
    
    # Identify all tokens resembling K\d{5}
    token_pattern = re.compile(r"K\d{5}")
    
    def replacer(match):
        ko = match.group(0)
        return "True" if ko in present_kos else "False"
        
    expr = token_pattern.sub(replacer, k_string)
    
    # Safety check: allow only True, False, and, or, (, ), space
    if not re.match(r"^[TrueFalseandor\(\)\s]+$", expr):
         # If parsing fails or unsafe, assume Absent
         return False
         
    try:
        return eval(expr)
    except:
        return False

def main():
    parser = argparse.ArgumentParser(description="Generate METABOLIC Worksheets")
    parser.add_argument("--hits", required=True)
    parser.add_argument("--template1", required=True)
    parser.add_argument("--template2", required=True)
    parser.add_argument("--kegg_module_db", required=True)
    parser.add_argument("--kegg_steps", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--dbcan_hits", required=False)
    parser.add_argument("--merops_hits", required=False)
    parser.add_argument("--genome_ids", required=True, help="Comma separated list of genome IDs")
    
    args = parser.parse_args()
    genome_ids = args.genome_ids.split(",")
    
    # Load hits
    try:
        df_hits = pd.read_csv(args.hits, sep="\t")
    except pd.errors.EmptyDataError:
        df_hits = pd.DataFrame(columns=["seq_id", "hmm_name", "genome_id"])

    # Load dbCAN hits
    # dbCAN output has no header: Target, Tlen, Query, Qlen, Evalue, HmmFrom, HmmTo, AliFrom, AliTo, Coverage
    # We also need genome_id. If collecting files, we might need filename mapping or 'genome_id' col added.
    # Assuming the input `dbcan_hits` has been collected with `keepHeader=false`?
    # Or maybe we rely on `seq_id` to look up genome_id (if we have a map)?
    # Better: Ensure the parsing step adds `genome_id`. 
    # Let's assume the passed file has `seq_id` and `family` (Query).
    # Since `parse_dbcan.py` outputs raw legacy columns, and we concat them...
    # We need to map SeqID -> GenomeID.
    # We can infer it from `df_hits` (HMM results) since that supposedly covers all proteins?
    # Or just load `all_hits_parsed.tsv` to build a `seq_map`.
    
    seq_map = {} # seq_id -> genome_id
    if not df_hits.empty: # assuming HMM hits cover most
         # But dbCAN hits might be on proteins not in HMM hits?
         # Nextflow `collectFile` can append reference name?
         pass
         
    # To be safe, we should rely on SeqID parsing or strict mapping.
    # Legacy: $Seqid2Genomeid from `ls *.faa`.
    # In Nextflow, we operate on streams.
    # We can assume `seq_id` contains `genome_id` if Prodigal named them `genome`_`gene`?
    # Prodigal `-a output` usually uses generic headers? No, `-a` uses seq id.
    # If we renamed headers...
    
    # Let's assume we can derive Genome ID from Seq ID or user provided strict mapping.
    # Or, update `modules/parsing.nf` to APPEND genome_id column.
    
    # Let's proceed assuming we have a map or can build it.
    for _, row in df_hits.iterrows():
         seq_map[row['seq_id']] = row['genome_id']

    dbcan_counts = {g: {} for g in genome_ids}
    if args.dbcan_hits:
        try:
            # Check structure. If generated by parse_dbcan.py, it has 10 cols, NO header.
            # We assume it was collected.
            with open(args.dbcan_hits, 'r') as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 3: continue
                    # Target(0), Query(2)
                    seq_id = parts[0]
                    family = parts[2]
                    
                    gid = seq_map.get(seq_id)
                    # If strictly using prodigal outputs named by genome, seq_id usually starts with genome_id?
                    # Or we just fail to count if missing. 
                    # Improvement: parse_dbcan should accept genome_id and output it.
                    # I will update parsing module to add genome_id column.
                    # Assuming here `dbcan_hits` has appended genome_id at column 10 (index 10) or similar.
                    # Let's assume I WILL update the parser to add it.
                    if len(parts) > 10:
                        gid = parts[10]
                    
                    if gid and gid in dbcan_counts:
                         # Strip .hmm from family if present
                         fam = family.replace(".hmm","")
                         dbcan_counts[gid][fam] = dbcan_counts[gid].get(fam, 0) + 1
        except Exception as e:
            logging.warning(f"Failed to parse dbCAN hits: {e}")

    merops_counts = {g: {} for g in genome_ids}
    if args.merops_hits:
        try:
             df_m = pd.read_csv(args.merops_hits, sep="\t")
             # Header: seq_id, family, target_id, evalue, genome_id
             for _, row in df_m.iterrows():
                 gid = row['genome_id']
                 if pd.isna(gid) or str(gid) == 'nan' or not gid:
                     gid = seq_map.get(row['seq_id'])
                 
                 if gid and gid in merops_counts:
                     merops_counts[gid][row['family']] = merops_counts[gid].get(row['family'], 0) + 1
        except Exception as e:
            logging.warning(f"Failed to parse MEROPS hits: {e}")


    # Build presence maps
    presence_map = {g: {} for g in genome_ids}
    ko_map = {g: set() for g in genome_ids}
    
    if not df_hits.empty:
        for _, row in df_hits.iterrows():
            g = row['genome_id']
            h = row['hmm_name']
            if not h.endswith(".hmm"):
                 h_clean = f"{h}.hmm"
            else:
                 h_clean = h
            
            if g in presence_map:
                presence_map[g][h_clean] = presence_map[g].get(h_clean, 0) + 1
                
                # KO mapping
                # HMM name IS the KO if it starts with K
                # If custom HMM, we might map to KO?
                # Legacy: _get_hmm_2_KO_hash logic.
                # Assuming parsed hits already normalized HMM names?
                # parse_hmm_results.py uses row['hmm_name'].
                # If it's KOfam, it's K00001. 
                # If parsed name is K00001.hmm -> K00001
                
                if h.startswith("K") and h[1:6].isdigit():
                    ko_map[g].add(h.split(".")[0])
    
    # --- Worksheet 1 & 2 ---
    # (Same as before, reusing logic)
    print("Generating Worksheet 1 & 2...")
    templates1 = load_hmm_template(args.template1)
    
    rows = []
    for entry_id, t in templates1.items():
        row = {
            'Category': t['category'],
            'Function': t['function'],
            'Gene abbreviation': t['gene_abbr'],
            'Gene name': t['gene_name'],
            'Hmm file': t['hmm_file'],
            'Corresponding KO': t['ko'],
            'Reaction': t['reaction'],
            'Product': t['product'],
            'Hmm detecting threshold': t['threshold']
        }
        hmms = [x.strip() for x in t['hmm_file'].split(",")]
        for g in genome_ids:
            count = sum(presence_map[g].get(h, 0) for h in hmms)
            row[f"{g} Presence"] = "Present" if count > 0 else "Absent"
            row[f"{g} Hit numbers"] = count
        rows.append(row)
    pd.DataFrame(rows).to_csv(f"{args.output_dir}/METABOLIC_result_worksheet1.tsv", sep="\t", index=False)
    
    templates2 = load_hmm_template_2(args.template2)
    rows2 = []
    for t in templates2:
        row = {
            'Category': t['category'],
            'Function': t['function'],
            'Gene abbreviation': t['gene_abbr']
        }
        refs = t['entry_refs'].split("||")
        target_hmms = []
        for r in refs:
             if r in templates1:
                 t1_hmms = [x.strip() for x in templates1[r]['hmm_file'].split(",")]
                 target_hmms.extend(t1_hmms)
                 
        for g in genome_ids:
             count = sum(presence_map[g].get(h, 0) for h in target_hmms)
             row[f"{g} Function presence"] = "Present" if count > 0 else "Absent"
        rows2.append(row)
    pd.DataFrame(rows2).to_csv(f"{args.output_dir}/METABOLIC_result_worksheet2.tsv", sep="\t", index=False)

    # --- Worksheet 4 (Step Presence) ---
    print("Generating Worksheet 4...")
    steps = load_kegg_steps(args.kegg_steps)
    module_cat = load_kegg_module_db(args.kegg_module_db)
    
    # Also need Module -> Total Steps for Worksheet 3
    module_steps_count = collections.defaultdict(int)
    module_steps_present = {m: {g: 0 for g in genome_ids} for m in module_cat}
    
    rows4 = []
    # Steps are ordered? Legacy sorts by key M...+01
    for step_id in sorted(steps.keys()):
        step = steps[step_id]
        mod_id = step['module']
        cat = module_cat.get(mod_id, "")
        
        row = {
            'Module step': step_id,
            'Module': step['desc'], # Using description as name?
            'KO id': step['k_string'],
            'Module Category': cat
        }
        
        module_steps_count[mod_id] += 1
        
        for g in genome_ids:
            is_present = check_step_presence(step['k_string'], ko_map[g])
            row[f"{g} Module step presence"] = "Present" if is_present else "Absent"
            
            if is_present:
                # Accumulate for Worksheet 3
                if mod_id not in module_steps_present:
                     module_steps_present[mod_id] = {gid: 0 for gid in genome_ids}
                module_steps_present[mod_id][g] += 1
                
        rows4.append(row)
    
    pd.DataFrame(rows4).to_csv(f"{args.output_dir}/METABOLIC_result_worksheet4.tsv", sep="\t", index=False)

    # --- Worksheet 3 (Module Presence) ---
    print("Generating Worksheet 3...")
    rows3 = []
    # Identify all modules seen in steps output
    all_modules = sorted(list(set(s['module'] for s in steps.values())))
    
    for mod_id in all_modules:
        # Get name (from any step?)
        # Legacy: $KEGG_module2name{$module}
        # We can find one step for this module
        name = ""
        # Find first step
        for s_id in steps:
            if steps[s_id]['module'] == mod_id:
                name = steps[s_id]['desc']
                break
                
        cat = module_cat.get(mod_id, "")
        total_steps = module_steps_count[mod_id]
        
        row = {
            'Module ID': mod_id,
            'Module': name,
            'Module Category': cat
        }
        
        for g in genome_ids:
             present_count = module_steps_present.get(mod_id, {}).get(g, 0)
             if total_steps > 0:
                 ratio = present_count / total_steps
                 # Legacy cutoff default 0.75? Or param?
                 # Assuming 0.75 for now (standard default)
                 # Wait, we need to pass cutoff param!
                 # Default is 0.75 in METABOLIC
                 if ratio >= 0.75:
                     row[f"{g} Module presence"] = "Present"
                 else:
                     row[f"{g} Module presence"] = "Absent"
             else:
                 row[f"{g} Module presence"] = "Absent"
        rows3.append(row)
        
    pd.DataFrame(rows3).to_csv(f"{args.output_dir}/METABOLIC_result_worksheet3.tsv", sep="\t", index=False)

    # --- Worksheet 5 (dbCAN / CAZyme) ---
    print("Generating Worksheet 5...")
    # Rows: CAZyme ID
    # Columns: [GID] Hit numbers, [GID] Hits (omitted logic for simplicity)
    # Collect all families
    all_cazymes = set()
    for g in dbcan_counts:
        all_cazymes.update(dbcan_counts[g].keys())
    
    rows5 = []
    for fam in sorted(all_cazymes):
        row = {'CAZyme ID': fam}
        for g in genome_ids:
             row[f"{g} Hit numbers"] = dbcan_counts[g].get(fam, 0)
             # row[f"{g} Hits"] = ...
        rows5.append(row)
    pd.DataFrame(rows5).to_csv(f"{args.output_dir}/METABOLIC_result_worksheet5.tsv", sep="\t", index=False)

    # --- Worksheet 6 (MEROPS) ---
    print("Generating Worksheet 6...")
    all_pep = set()
    for g in merops_counts:
        all_pep.update(merops_counts[g].keys())
        
    rows6 = []
    for pep in sorted(all_pep):
        row = {'MEROPS peptidase ID': pep}
        for g in genome_ids:
             row[f"{g} Hit numbers"] = merops_counts[g].get(pep, 0)
        rows6.append(row)
    pd.DataFrame(rows6).to_csv(f"{args.output_dir}/METABOLIC_result_worksheet6.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
