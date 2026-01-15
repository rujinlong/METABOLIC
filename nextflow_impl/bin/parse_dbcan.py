#!/usr/bin/env python3

import sys
import argparse
import logging
import csv

# Legacy defaults from hmmscan-parser-dbCANmeta.py
DEFAULT_EVAL = 1e-15
DEFAULT_COVERAGE = 0.35

def parse_domtblout_line(line):
    """
    Parse a domtblout line (whitespace delimited).
    Standard Domtblout cols:
    0: target name (sequence)
    1: accession
    2: tlen
    3: query name (hmm)
    4: accession
    5: qlen
    6: E-value
    ...
    15: from (hmm)
    16: to (hmm)
    17: from (ali/seq)
    18: to (ali/seq)
    19: from (env)
    20: to (env)
    """
    parts = line.strip().split()
    if len(parts) < 22: return None
    
    # We need:
    # HMM Name (Query) -> col 3
    # Seq ID (Target) -> col 0
    # E-value -> col 6
    # HMM Len -> col 5
    # Seq From -> col 17
    # Seq To -> col 18
    # Seq Len -> col 2 (tlen)
    
    # Also need HMM start/end for coverage? 
    # Legacy script uses:
    # awk $1,$3,$4,$6,$13,$16,$17,$18,$19
    # $1=Target Name
    # $3=Target Len (tlen)
    # $4=Query Name (HMM)
    # $6=Query Len (qlen)
    # $13=i-Evalue (or E-value? Usually col 6 is full seq evalue, 12 is i-Evalue? wait)
    # Let's check HMMER manual or assume legacy awk indices are correct 1-based.
    # Awk $1 = Col 0
    # Awk $3 = Col 2
    # Awk $4 = Col 3
    # Awk $6 = Col 5
    # Awk $13 = Col 12 (i-Evalue domain)
    # Awk $16 = Col 15 (hmm from)
    # Awk $17 = Col 16 (hmm to)
    # Awk $18 = Col 17 (ali from)
    # Awk $19 = Col 18 (ali to)
    
    return {
        'target_name': parts[0],
        'tlen': int(parts[2]),
        'hmm_name': parts[3],
        'qlen': int(parts[5]),
        'evalue': float(parts[12]), # Using i-Evalue as per legacy index 13
        'hmm_from': int(parts[15]),
        'hmm_to': int(parts[16]),
        'ali_from': int(parts[17]),
        'ali_to': int(parts[18])
    }

def resolve_overlaps(hits):
    """
    Resolve overlapping hits for a single QUERY(HMM)?
    No, legacy groups by "$a[2]" -> Column 3 -> Query Name -> HMM.
    Wait...
    Legacy script: `awk '{print $1,...}' | ... push(@{$b{$a[2]}},$_)`
    Awk output: 
    0: $1 Target(Seq)
    1: $3 TLen
    2: $4 Query(HMM)
    
    So legacy groups by HMM Name.
    Then it sorts by Target(Seq) (sort -k 3,3 ?? Wait. 
    Awk output index 2 is Query(HMM).
    Awk output index 0 is Target(Seq).
    
    Legacy Sort: `-k 3,3` -> Column 3?
    The output of awk is space delimited (sed replaces space with tab).
    Columns sent to perl:
    0: Target
    1: Tlen
    2: Query(HMM)
    ...
    
    So it sorts by Query(HMM).
    Then Perl code: `push(@{$b{$a[2]}},$_)` -> Keys are Query(HMM).
    Wait, `sort -k 3,3` sorts by 3rd column designated by space/tab?
    If input to sort has tabs... yes.
    
    So it processes hits PER HMM.
    And within HMM, it checks overlaps?
    `@a=@{$b{$_}}` (Hits for one HMM).
    `for ($i=0...)` compare $a[$i] and $a[$i+1].
    The hits are sorted by Col 8 and 9 (Ali From/To) `-k 8n -k 9n`.
    (Awk output index 7 and 8).
    
    So: For a given HMM, sort hits by Start coordinate on Sequence?
    Wait. `hmmscan` runs on ONE sequence vs DB? Or DB vs Seq?
    HMMScan: Search sequence(s) against HMM DB.
    So different sequences are independent.
    If we group by HMM, we are mixing hits from different sequences?
    Unless the sort includes sequence ID?
    Legacy sort: `-k 3,3` (HMM Name).
    It does NOT sort by Sequence ID. 
    So it mixes sequences if multiple sequences hit the same HMM.
    
    BUT `next if $a[-1]==$a[-2]` (Skip if length 0?)
    
    Then: `if ($len3>0 ...)` -> overlap calculation.
    `$b` and `$c` are adjacent hits in the sorted list.
    If they are different sequences, their coordinates are incomparable.
    Does legacy buggily filter across sequences?
    Or maybe `hmmscan` output is grouped by sequence already?
    
    Standard `hmmscan` output groups by target sequence if not sorted.
    But legacy creates a hash `%b` keyed by HMM Name.
    Then it iterates values in that hash.
    So it compares hits for HMM X on Seq A with hits for HMM X on Seq B?
    Coordinates might overlap by chance.
    
    This looks like a BUG in legacy code or I misunderstand the column indices.
    Awk: $1(Target), $3(Tlen), $4(Query)...
    Col 2 in Perl (0-based) is output index 2 -> Query.
    So yes, keyed by HMM.
    
    If Target A has hit 100-200.
    Target B has hit 100-200.
    Are they considered overlapping?
    If they are sorted by coords only...
    
    Wait, `sort -k 3,3`. This sorts by Query.
    So all hits for same HMM are together.
    Inside Perl: `foreach(sort keys %b)` -> Iterate HMMs.
    `@a` = list of hits for this HMM.
    But `@a` is populated from the sorted stream. The sort didn't sort by Sequence ID within HMM?
    Actually `sort -k 3,3 -k 8n -k 9n`.
    Primary key: HMM. Secondary: Start. Tertiary: End.
    It does NOT sort by Sequence.
    
    So if Seq A (Hit 100-200) and Seq B (Hit 150-250) are both in the list.
    Adjacent? Yes.
    Overlap? Yes.
    Legacy code removes one?
    `if ($b[4] < $c[4])` (Compare E-value/Score? Col 4 in output is $13 in awk -> i-Evalue).
    So it keeps better E-value.
    
    This implies Legacy code effectively allows only ONE hit per HMM across the entire dataset if they happen to have overlapping LOCAL coordinates?
    Or maybe `dbCAN` assumes run on single genome?
    Legacy loop: `while (<IN_FAA>) ... hmmscan ... $gn_id.dbCAN2.out`.
    It runs PER GENOME text file.
    Inside that file, seq IDs are different.
    So `genome_id` is constant. But `seq_id` (Gene) varies.
    
    If Gene A and Gene B both have hits to HMM X at same coords (e.g. 10-100), they will be seen as overlapping and one removed?
    This is extremely aggressive/buggy if so.
    However, I must assume "Result Consistency" is paramount.
    But "Fixed Bug" is allowed if "Unnecessary intermediate files" rule applies? No.
    "final output results to be consistent".
    
    Let's check if `hmmscan` output `domtblout` usually includes multiple genes per file. Yes.
    
    Hypothesis: The Legacy `hmmscan-parser` is intended to filter overlapping domains ON THE SAME SEQUENCE.
    But the implementation groups by HMM and ignores Sequence ID during overlap check.
    This means if Gene A (1-100) and Gene B (1-100) both hit HMM X, they are treated as overlapping and one is deleted.
    Since gene coordinates usually start at 1 for Profigal... almost ALL genes overlap in local coords.
    So this effectively keeps ONLY ONE HIT PER HMM per Genome (mostly).
    
    Is this desired?
    dbCAN usually counts CAZymes.
    If I have 5 copies of GH1, I want to count 5.
    If this parser removes 4 because they all map to 1-300 aa...
    Then the result is wrong.
    
    However, maybe `$a[2]` is NOT HMM name?
    Awk `$1,$3,$4`.
    $1 = Target (Seq). $4 = Query (HMM).
    Perl `split; ... $a[2]`.
    Perl split on whitespace.
    The awk output is space separated (then `sed` to tab? `sed 's/ /\t/g'`).
    So Perl sees tab separated.
    Index 0: Target
    Index 1: Tlen
    Index 2: Query (HMM)
    
    Yes, it keys by HMM.
    
    I will replicate this logic exactly (even if it seems buggy) to ensure consistency, 
    UNLESS it's obviously just wrong interpretation of my part.
    
    Wait, `hmmscan-parser-dbCANmeta.py` is called on the output of `hmmscan`.
    The inputs to `hmmscan` in `METABOLIC-G.pl` are `$file` (individual `.faa` per genome).
    So result is per-genome.
    
    If I run this in Python, I should probably group by Target(Seq) AND Query(HMM)?
    No, usually domain filtering is Per-Sequence, Per-HMM (or Per-Sequence, All-Domains).
    
    If I group by HMM, I am filtering domains for that HMM across all sequences.
    
    Let's check `sort -k 3,3`.
    If I assume `$a[2]` is actually Target Name?
    If awk printed `$4` then `$1`...
    But awk prints `$1,$3,$4`.
    
    Maybe I should verify what `$1` and `$4` are in `hmmscan --domtblout`.
    Reference: HMMER User Guide.
    domtblout:
    1. target name
    2. accession
    3. tlen
    4. query name
    
    So $1=target, $4=query.
    
    Conclusion: Legacy script groups by Query (HMM).
    So on a per-genome basis, it seems to filter hits to the same HMM that share coordinates, regardless of which gene they are on.
    Since genes start at 0/1, they share coordinates.
    So it essentially picks the BEST gene for each HMM if they have similar lengths?
    
    This effectively "Deduplicates" HMM hits so each CAZyme family appears only once (or few times) per genome?
    If true, Worksheet 5 should show very low counts.
    
    Let's replicate it faithfully.
    
    Algorithm:
    1. Parse all hits.
    2. Group by HMM Name.
    3. Sort by Ali_From, Ali_To.
    4. Iterate and remove overlaps (using the len3 logic).
    5. Filter by coverage and evalue.
    
    Wait, `row.append(float(int(row[6])-int(row[5]))/int(row[1]))`
    `row[6]` is Index 6 in Python list (which matches Index 6 in Perl/Awk stream?).
    Lines coming to Python `with open('temp_')`.
    The perl script prints filtered/spliced array `@a`.
    It prints `join`ed string? `print $_`.
    It preserves input format (Tab separated 9 items).
    0: Target
    1: Tlen
    2: Query
    3: Qlen
    4: Evalue
    5: HmmFrom
    6: HmmTo
    7: AliFrom
    8: AliTo
    
    Python calculation: `(row[6] - row[5]) / row[1]`.
    (HmmTo - HmmFrom) / Tlen(SeqLen).
    This is "Coverage of Sequence"? No.
    (Hmm coordinates) / (Sequence Length).
    Usually coverage is (AliTo - AliFrom) / (HmmLen) for "HMM Coverage".
    Or (AliTo - AliFrom) / (SeqLen) for "Seq Coverage".
    
    Here it takes HMM coords (on the model) divided by Target Sequence Length?
    That seems... odd.
    Usually you divide by Qlen (HMM Length) to see how much of model is covered.
    Or divide by Tlen to see how much of gene is domain.
    
    But `(HmmTo - HmmFrom)` is the length covered on the MODEL.
    `row[1]` is `Tlen` (Target Sequence Length).
    So it calculates "Fraction of Gene covered by Model State Path".
    
    Wait, let's verify indices.
    Awk: `$1,$3,$4,$6,$13,$16,$17,$18,$19`
    0: Target
    1: Tlen
    2: Query
    3: Qlen
    4: Evalue
    5: HmmFrom (16)
    6: HmmTo (17)
    7: AliFrom (18)
    8: AliTo (19)
    
    Python `row[6]` is `HmmTo`. `row[5]` is `HmmFrom`. `row[1]` is `Tlen`.
    Yes. `(HmmTo - HmmFrom) / Tlen`.
    
    And check `coverage >= 0.35`.
    
    Okay, I will document this oddity but implement EXACTLY this.

    Usage: `python parse_dbcan.py input.dm output.tsv`
    """
    pass  # End of docstring/comment function

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="domtblout file")
    parser.add_argument("output", help="output file")
    parser.add_argument("--genome_id", required=True, help="Genome ID to append to output")
    args = parser.parse_args()
    
    hits_by_hmm = {}
    
    # 1. Parse and Group
    with open(args.input) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 22: continue
            
            # Extract necessary fields mapping to legacy
            # target=0, tlen=2, query=3, qlen=5, eval=12, hmmfrom=15, hmmto=16, alifrom=17, alito=18
            
            item = {
                'line': line, # Keep original line? No, legacy outputs specific cols
                'target': parts[0],
                'tlen': int(parts[2]),
                'query': parts[3],
                'qlen': int(parts[5]),
                'evalue': float(parts[12]),
                'hmm_from': int(parts[15]),
                'hmm_to': int(parts[16]),
                'ali_from': int(parts[17]),
                'ali_to': int(parts[18])
            }
            
            if item['query'] not in hits_by_hmm:
                hits_by_hmm[item['query']] = []
            hits_by_hmm[item['query']].append(item)
            
    # 2. Filter Overlaps
    final_hits = []
    
    for hmm, hits in hits_by_hmm.items():
        # Sort by AliFrom (asc), AliTo (asc)
        # Perl: sort -k 8n -k 9n (indices 7,8 0-based from awk) -> AliFrom, AliTo
        hits.sort(key=lambda x: (x['ali_from'], x['ali_to']))
        
        # Iterative filtering
        i = 0
        while i < len(hits) - 1:
            curr = hits[i]
            next_h = hits[i+1]
            
            # Lengths on sequence
            len1 = curr['ali_to'] - curr['ali_from']
            len2 = next_h['ali_to'] - next_h['ali_from']
            
            # Overlap
            # Since sorted by from: overlap is min(to) - next_from
            overlap_end = min(curr['ali_to'], next_h['ali_to'])
            overlap_start = next_h['ali_from']
            len3 = overlap_end - overlap_start
            
            # Perl logic: $len3 > 0 and ($len3/$len1 > 0.5 or $len3/$len2 > 0.5)
            if len3 > 0:
                ratio1 = len3 / len1 if len1 > 0 else 0
                ratio2 = len3 / len2 if len2 > 0 else 0
                
                if ratio1 > 0.5 or ratio2 > 0.5:
                    # Remove the one with worse E-value
                    # Perl: if ($b[4] < $c[4]) -> b is smaller (better) -> remove c
                    # b=curr (i), c=next (i+1). 4=Evalue.
                    if curr['evalue'] < next_h['evalue']:
                        # Remove next
                        hits.pop(i+1)
                        # Don't increment i, check new next
                    else:
                        # Remove curr
                        hits.pop(i)
                        # i stays same (which is now new curr)
                        # But we need to handle index underflow?
                        # `i` points to current slot. we popped it. now `i` points to old `next`.
                        # We should re-check this new item against ITS next.
                        # BUT Perl logic: `$i = $i - 1`.
                        if i > 0:
                            i -= 1
                        else:
                            # if i=0, we stay at 0.
                            pass
                    continue
            i += 1
            
        final_hits.extend(hits)
        
    # 3. Filter by Coverage and Write
    with open(args.output, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        # Legacy output columns: 
        # Target, Tlen, Query, Qlen, Evalue, HmmFrom, HmmTo, AliFrom, AliTo, Coverage
        
        for h in final_hits:
            # Coverage calculation: (HmmTo - HmmFrom) / Tlen
            # float(int(row[6])-int(row[5]))/int(row[1])
            if h['tlen'] == 0: cov = 0
            else: cov = (h['hmm_to'] - h['hmm_from']) / h['tlen']
            
            if h['evalue'] <= DEFAULT_EVAL and cov >= DEFAULT_COVERAGE:
                writer.writerow([
                    h['target'], h['tlen'], h['query'], h['qlen'], h['evalue'],
                    h['hmm_from'], h['hmm_to'], h['ali_from'], h['ali_to'], cov,
                    args.genome_id
                ])

if __name__ == "__main__":
    main()
