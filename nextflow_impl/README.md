# METABOLIC-Nextflow Implementation

A streamlined, HPC-optimized Nextflow implementation of the **METABOLIC** (Metabolic And Biogeochemistry Analyses In Microbes) pipeline.

## 🚀 Key Features

- **Single Search Architecture**: One `hmmsearch` call for all proteins (vs. N per genome)
- **MMseqs2 Clustering**: Optional protein dereplication reduces search time by 50-80%
- **HPC-Optimized**: Designed for large memory (500GB+) and high CPU (50+) environments
- **Pre-merged HMM Databases**: Eliminates per-HMM I/O overhead
- **Post-Search Thresholding**: Python-based filtering with per-HMM thresholds

---

## 📥 Input Requirements

### Required Files

1. **`proteins.faa`** - Concatenated protein sequences from all genomes/MAGs
2. **`contig_map.tsv`** - Maps contig IDs to genome/MAG IDs

### Protein FASTA Format

> [!IMPORTANT]
> Protein headers **MUST** follow this format: `>contigID_proteinID`

```fasta
>contig001_1 hypothetical protein
MKKLLLLLLLLLLLLLL...
>contig001_2 ATP synthase
MVVVVVVVVVVVVVVVV...
>contig002_1 ribosomal protein
MAAAAAAAAAAAAAAAA...
```

**Naming Convention:**
- `contigID`: The contig identifier (must match `contig_map.tsv`)
- `proteinID`: Integer suffix for protein number on that contig
- Separator: Underscore (`_`) between contigID and proteinID

**If using Prodigal/Pyrodigal:**
```bash
# Prodigal default output already uses this format
prodigal -i genome.fna -a proteins.faa
# Output: >contig001_1, >contig001_2, etc.
```

### Contig Mapping File Format

Tab-separated, no header:
```tsv
contig001	MAG_001
contig002	MAG_001
contig003	MAG_002
scaffold_1	Genome_A
scaffold_2	Genome_A
```

**Column 1**: Contig ID (must match prefix in protein headers)
**Column 2**: Genome/MAG ID

---

## 💻 Usage

### Basic Run (with Clustering)
```bash
nextflow run main.nf \
    --proteins /path/to/all_proteins.faa \
    --contig_map /path/to/contig_to_genome.tsv \
    --outdir ./results \
    -profile singularity
```

### Skip Clustering (Pre-dereplicated Input)
```bash
nextflow run main.nf \
    --proteins /path/to/dereplicated_proteins.faa \
    --contig_map /path/to/contig_to_genome.tsv \
    --skip_cluster \
    -profile singularity
```

### Custom Clustering Parameters
```bash
nextflow run main.nf \
    --proteins /path/to/all_proteins.faa \
    --contig_map /path/to/contig_to_genome.tsv \
    --cluster_sid 0.95 \
    --cluster_cov 0.9 \
    -profile hpc
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--proteins` | (required) | Path to concatenated protein FASTA |
| `--contig_map` | (required) | Path to contig-genome mapping TSV |
| `--outdir` | `results` | Output directory |
| `--db_dir` | `$HOME/database/METABOLIC` | Database root directory |
| **Clustering** | | |
| `--skip_cluster` | `false` | Skip MMseqs2 clustering |
| `--cluster_sid` | `0.9` | Minimum sequence identity (90%) |
| `--cluster_cov` | `0.8` | Minimum coverage (80%) |

---

## 🎯 Custom DB Priority

By default, when the same protein hits both **custom_db** and **KOfam** for the same underlying function, the pipeline **prefers custom_db hits**.

### Why?

- **Custom HMMs are more specific**: Curated for biogeochemical cycles in METABOLIC
- **Avoid double-counting**: Same gene shouldn't be counted twice
- **Better motif validation**: Custom HMMs often have stricter active-site requirements

### How It Works

1. The `hmm_table_template.txt` contains mappings: `custom_hmm` → `corresponding_KO`
2. When a protein matches both `amoA.hmm` (custom) and `K10944` (KOfam):
   - Both hits pass threshold filtering
   - Deduplication removes the KOfam hit (since custom covers K10944)
   - Only the custom hit is kept

### Behavior Control

The `--prefer_custom` flag is **enabled by default**. To disable:

```bash
# Keep all hits (no deduplication)
nextflow run main.nf ... --no-prefer-custom
```

### Example Output

```
[12:00:03] Parsed 150000 raw hits (KOfam + Custom)
[12:00:05] Filtered to 45000 valid hits
[12:00:05] Custom DB priority: removed 823 redundant KOfam hits
[12:00:05] After deduplication: 44177 hits
```

---

## 💾 Database Setup

### Directory Structure
```
${db_dir}/
├── kofam_database/
│   ├── profiles/           # Original HMMs (optional, for reference)
│   ├── kofam_all.hmm       # Merged KOfam HMMs
│   └── ko_list             # Thresholds from KEGG
├── custom_db/
│   └── metabolic_all.hmm   # Merged custom METABOLIC HMMs
├── METABOLIC_template_and_database/
│   ├── hmm_table_template.txt
│   ├── hmm_table_template_2.txt
│   ├── ko00002.keg
│   ├── kegg_module_step_db.txt
│   ├── motif.txt
│   └── motif.pair.txt
├── dbCAN2/
│   └── dbCAN-fam-HMMs.txt  # Must be indexed with hmmpress
└── MEROPS/
    ├── pepunit.db.dmnd     # Diamond database
    └── pepunit.lib
```

### Merge Commands (One-time Setup)

```bash
DB_DIR="${HOME}/database/METABOLIC"

# 1. Merge KOfam HMMs (~21,000 files)
cat ${DB_DIR}/kofam_database/profiles/*.hmm > ${DB_DIR}/kofam_database/kofam_all.hmm

# 2. Merge Custom METABOLIC HMMs (normalize versions, fix NAME fields)
mkdir -p ${DB_DIR}/custom_db

for f in ${DB_DIR}/METABOLIC_hmm_db/*.hmm; do 
  name=$(basename "$f" .hmm)
  hmmconvert "$f" | sed "s/^NAME  .*/NAME  $name/"
done > ${DB_DIR}/custom_db/metabolic_all.hmm

# 3. Index merged HMMs (required for hmmscan/dbCAN)
hmmpress ${DB_DIR}/custom_db/metabolic_all.hmm
hmmpress ${DB_DIR}/dbCAN2/dbCAN-fam-HMMs.txt
```

---

## 🔬 Pipeline Workflow

```mermaid
graph TD
    A[proteins.faa] --> B{skip_cluster?}
    B -->|No| C[MMSEQS_CLUSTER]
    B -->|Yes| D[proteins]
    C --> D[rep_seq.faa]
    C --> E[cluster.tsv]
    
    D --> F[HMMSEARCH_KO]
    D --> G[HMMSEARCH_CUSTOM]
    D --> H[DBCAN_SEARCH]
    D --> I[MEROPS_SEARCH]
    
    F --> J[PARSE_HMM]
    G --> J
    E --> J
    H --> K[PARSE_DBCAN]
    E --> K
    I --> L[PARSE_MEROPS]
    E --> L
    
    J --> M[GENERATE_TABLES]
    K --> M
    L --> M
    
    M --> N[Worksheet 1-6]
```

### Steps

1. **Cluster** (optional): MMseqs2 clusters proteins by 90% identity (default)
2. **Search**: Representative proteins searched against merged HMM databases
3. **Parse**: Filter hits by thresholds, expand clusters, resolve genome IDs
4. **Generate**: Produce 6 worksheets summarizing metabolic potential

---

## 📊 Output Files

All output files are TSV format with snake_case column names.

### Intermediate Files

| File | Description |
|------|-------------|
| `all_hits_parsed.tsv` | All filtered HMM hits with genome assignments |
| `all_dbcan_parsed.tsv` | CAZyme domain hits |
| `all_merops_parsed.tsv` | Peptidase hits |

#### `all_hits_parsed.tsv` Columns

| Column | Description |
|--------|-------------|
| `seq_id` | Protein ID (e.g., `BBKHAN_02371`) |
| `hmm_name` | HMM that matched (e.g., `K00003` or `amoA`) |
| `full_score` | Full sequence bit score |
| `domain_score` | Best domain bit score |
| `source` | Hit source: `kofam` or `custom` |
| `genome_id` | Genome/MAG ID from contig_map |

#### `all_dbcan_parsed.tsv` Columns

| Column | Description |
|--------|-------------|
| `seq_id` | Protein ID |
| `hmm_name` | CAZyme family (e.g., `GH1`, `AA3`) |
| `hmm_len` | HMM length |
| `seq_len` | Protein length |
| `evalue` | Domain i-Evalue |
| `coverage` | HMM coverage |
| `genome_id` | Genome/MAG ID |

---

### Final Worksheets

#### Worksheet 1: Gene-level Hits

Detailed gene presence/absence with hit counts per genome.

| Column | Description |
|--------|-------------|
| `category` | Functional category |
| `function` | Metabolic function |
| `gene_abbr` | Gene abbreviation |
| `gene_name` | Full gene name |
| `hmm_file` | HMM file(s) used |
| `corresponding_ko` | KEGG Orthology ID |
| `reaction` | KEGG reaction |
| `product` | Reaction product |
| `hmm_threshold` | Detection threshold |
| `{genome}_presence` | Present/Absent |
| `{genome}_hit_count` | Number of hits |

#### Worksheet 2: Function Summary

Aggregated function presence per genome.

| Column | Description |
|--------|-------------|
| `category` | Functional category |
| `function` | Metabolic function |
| `gene_abbr` | Gene abbreviation |
| `{genome}_function_presence` | Present/Absent |

#### Worksheet 3: KEGG Module Completeness

Module-level presence (≥75% steps = Present).

| Column | Description |
|--------|-------------|
| `module_id` | KEGG Module ID (e.g., `M00001`) |
| `module_name` | Module name |
| `module_category` | Module category |
| `{genome}_module_presence` | Present/Absent |

#### Worksheet 4: KEGG Module Steps

Step-by-step module presence.

| Column | Description |
|--------|-------------|
| `module_step` | Step ID (e.g., `M00001+01`) |
| `module_name` | Module name |
| `ko_expression` | Boolean KO expression |
| `module_category` | Module category |
| `{genome}_step_presence` | Present/Absent |

#### Worksheet 5: CAZyme Summary

CAZyme family hit counts per genome.

| Column | Description |
|--------|-------------|
| `cazyme_id` | CAZyme family (e.g., `GH1`) |
| `{genome}_hit_count` | Number of hits |

#### Worksheet 6: MEROPS Summary

Peptidase family hit counts per genome.

| Column | Description |
|--------|-------------|
| `merops_id` | MEROPS family (e.g., `S01`) |
| `{genome}_hit_count` | Number of hits |

---

## ⚙️ Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Local execution |
| `docker` | Docker container |
| `singularity` | Singularity container |
| `hpc` | SLURM cluster (high memory/CPUs) |
| `test` | Test with sample data |
