# LemonIte PKN Build Pipeline

A three-step Python pipeline for constructing a **Prior Knowledge Network (PKN)** that integrates metabolite–gene and protein–protein interactions from nine curated databases. The PKN is used as a causal scaffold for multi-omics mechanistic modelling in the LemonIte project.

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Prerequisites](#prerequisites)
4. [Installation](#installation)
5. [Configuration](#configuration)
6. [Running the Pipeline](#running-the-pipeline)
7. [Database Sources](#database-sources)
8. [Output Files](#output-files)
9. [Caching System](#caching-system)
10. [HPC Deployment](#hpc-deployment)
11. [Singularity Container](#singularity-container)
12. [Troubleshooting](#troubleshooting)
13. [Repository Layout](#repository-layout)

---

## Overview

The pipeline builds a heterogeneous network with two edge types:

| Edge type | Meaning |
|---|---|
| `metabolite-gene` | A metabolite regulates / interacts with a gene product |
| `PPI` | Protein–protein interaction between two gene products |

Each edge carries a `Source` label (which database it was derived from), a `url` field pointing to the source record, and optionally a `reaction_id` (for metabolic model edges).

The final output, `LemonIte_PKN_with_URLs.tsv`, is ready for use as input to CORNETO, COSMOS, or any causal reasoning tool that accepts a signed interaction network.

---

## Architecture

```
PKN_build_pipeline/
├── main.py                  ← CLI entry point; orchestrates all steps
├── config.py                ← All paths, thresholds, and API settings
├── config_hpc.py            ← HPC-specific overrides (VSC Ghent)
│
├── step1_metabolites/       ← Metabolite–gene interaction retrieval
│   ├── preprocessing.py     ← HMDB loading, canonical SMILES, ChEMBL mapping
│   ├── biogrid.py           ← BioGRID chemical–protein interactions (local file)
│   ├── stitch.py            ← STITCH protein–chemical links (local file + MyGene API)
│   ├── uniprot.py           ← UniProtKB metabolite–protein (REST API)
│   ├── intact.py            ← IntAct metabolite–protein (REST API, ChEBI-based)
│   ├── chembl.py            ← ChEMBL bioactivity (REST API)
│   ├── lincs.py             ← LINCS biochemical binding (local files, IC50 threshold)
│   ├── gem.py               ← Human-GEM metabolic model (NetworkX graph traversal)
│   ├── metalinks.py         ← MetalinksDB ligand–receptor interactions (local file)
│   └── integration.py       ← Merges all step-1 results → metabolite_gene_PKN.tsv
│
├── step2_proteins/          ← Protein–protein interaction retrieval
│   ├── biogrid_ppi.py       ← BioGRID PPI (local file)
│   ├── huri.py              ← HuRI human interactome (local file)
│   ├── humannet.py          ← HumanNet PPI (local file)
│   ├── string_api.py        ← STRING v12 (REST API, confidence-filtered)
│   └── ppi_integration.py   ← Merges all step-2 results → PPI_network.tsv
│
├── step3_final/             ← Final PKN assembly
│   ├── combiner.py          ← Merges metabolite-gene + PPI networks
│   ├── annotator.py         ← Reports URL coverage; writes LemonIte_PKN_with_URLs.tsv
│   ├── analysis.py          ← Network statistics and coverage analysis
│   └── visualization.py     ← Plots (upset, barplot, degree distribution)
│
└── utils/
    ├── pipeline.py          ← Abstract base classes: DatabaseRetriever,
    │                           LocalFileRetriever, APIRetriever
    ├── api_retry.py         ← Exponential-backoff retry decorator (HTTP 429 aware)
    └── file_io.py           ← HMDB XML parser, ChEBI/PubChem extraction helpers
```

### Data flow

```
HMDB XML
   │
   ▼
preprocessing.py   ──── canonical SMILES, ChEMBL IDs ────┐
   │                                                       │
   ▼                                                       │
Step 1 Retrievers (biogrid, stitch, uniprot,               │
                   intact, chembl, lincs,                  │
                   gem_dist1, gem_dist2, metalinks)        │
   │                                                       │
   │  Each retriever writes:                               │
   │    HMDB_metabolites_<DB>_processed.csv  (cache)       │
   │    columns: [HMDB_ID, Gene, Source, url,              │
   │              reaction_id (GEM only)]                  │
   │                                                       │
   ▼                                                       │
integration.py  ──────────────────────────────────────────┘
   │  metabolite_gene_PKN.tsv
   │  columns: [Metabolite, Gene, Source, url, reaction_id]
   │
   ▼
Step 2 Retrievers (biogrid_ppi, huri, humannet, string)
   │  PPI_network.tsv
   │  columns: [GeneA, GeneB, Source, combined_score]
   │
   ▼
combiner.py
   │  LemonIte_PKN.tsv
   │  columns: [Node1, Node2, Type, Source, url, reaction_id, combined_score]
   │
   ▼
annotator.py  (reports URL coverage, passes url column through)
   │  LemonIte_PKN_with_URLs.tsv  ← final output
   ▼
analysis.py + visualization.py
```

---

## Prerequisites

### Required local database files

| Database | File | Used in |
|---|---|---|
| HMDB | `hmdb_metabolites.xml` | preprocessing (step 1) |
| BioGRID chemicals | `BIOGRID-CHEMICALS-*.chemtab.txt` | step 1 – BioGRID |
| STITCH | `9606.protein_chemical.links.v5.0.tsv` | step 1 – STITCH |
| Human-GEM | `Human-GEM.txt`, `metabolites.tsv`, `genes.tsv` | step 1 – GEM |
| MetalinksDB | `metalinks.csv` | step 1 – MetalinksDB |
| LINCS | `lsp_compound_mapping.csv`, `lsp_biochem_agg.csv`, `lsp_target_dictionary.csv` | step 1 – LINCS |
| BioGRID PPI | `BIOGRID-ALL-*.tab3.txt` | step 2 – BioGRID PPI |
| HuRI | `HuRI_Tong_2021.tsv` | step 2 – HuRI |
| HumanNet | `HS-PI.symbol.tsv` | step 2 – HumanNet |
| Ensembl mapping | `ensembl_to_hgnc.tsv` | step 2 – HuRI |

### External APIs (require internet access)

| API | Used for | Rate limiting |
|---|---|---|
| UniProtKB REST | Metabolite–protein interactions (InChIKey search) | max 3 threads, pause every 50 requests |
| IntAct EBI REST | Metabolite–protein interactions (ChEBI ID search) | max 4 threads, pause every 75 requests |
| ChEMBL REST | Bioactivity data | max 5 threads, pause every 100 requests |
| ChEMBL WebClient | SMILES → ChEMBL ID mapping (preprocessing) | max 5 threads |
| STRING v12 API | PPI network (gene-set query) | chunked requests |
| MyGene.info | Ensembl protein ID → gene symbol (STITCH) | batch queries, chunks of 1000 |

---

## Installation

### Local / interactive

```bash
cd PKN_build_pipeline
bash setup.sh          # creates venv, installs requirements.txt
source venv/bin/activate
```

### HPC (VSC Ghent)

A pre-built virtual environment is available at:
```
/user/gent/435/vsc43501/25BVM_Lemonite/venv_pkn_doduo
```

Activate it before running:
```bash
source /user/gent/435/vsc43501/25BVM_Lemonite/venv_pkn_doduo/bin/activate
```

### Dependencies

```
pandas>=2.0.0          numpy>=1.24.0
requests>=2.31.0       chembl-webresource-client>=0.10.8
networkx>=3.1          mygene>=3.2.2
matplotlib>=3.7.0      seaborn>=0.12.0
tqdm>=4.65.0           openpyxl>=3.1.0
lxml>=4.9.0            upsetplot>=0.8.0
matplotlib-venn>=0.11.0 rdkit>=2023.3.1
```

---

## Configuration

All settings live in `config.py`. The most common overrides are set via **environment variables** so that container / HPC runs do not require editing the file:

| Environment variable | Default | Description |
|---|---|---|
| `PKN_WORKDIR` | `/home/borisvdm/Documents/PhD/Lemonite` | Root working directory |
| `PKN_OUTPUT_DIR_NAME` | `PKN` | Sub-directory for all output files |
| `PKN_DB_DIR` | `/run/media/.../databases` | Root directory for all local database files |
| `PKN_GEM_DIR` | — | Directory containing Human-GEM model files |
| `PKN_CONFIG` | `config` | Python module name to use as config (e.g. `config_hpc`) |

HPC runs use `config_hpc.py` (loaded automatically by `submit_pkn_hpc.sh`), which sets correct HPC paths without modifying `config.py`.

### Key thresholds (in `config.py`)

| Parameter | Default | Meaning |
|---|---|---|
| `PATHWAY_DISTANCE` | `2` | GEM: max reaction hops from metabolite to gene |
| `GEM_UBIQUITOUS_METABOLITES` | list of 17 IDs | Currency metabolites excluded from graph traversal |
| STRING confidence | `400` | STRINGRetriever minimum combined score (medium confidence) |
| LINCS IC50 threshold | `10 000 nM` | LINCSRetriever maximum IC50 (10 µM) |

### API retry settings

Each API database has its own retry profile in `API_RETRY_CONFIG`:

```python
'UniProtKB': {
    'max_retries': 10,   # attempts before giving up
    'backoff_factor': 3, # wait = backoff_factor^(attempt+2) seconds
    'timeout': 30,       # per-request timeout (s)
    'max_workers': 3,    # concurrent threads
    'pause_after': 50,   # periodic pause every N requests
    'pause_duration': 15 # seconds to pause
}
```

---

## Running the Pipeline

### Full run

```bash
python main.py --all
```

### Individual steps

```bash
python main.py --step 1   # metabolite–gene interactions
python main.py --step 2   # protein–protein interactions
python main.py --step 3   # final PKN assembly + annotation
```

### Useful options

```bash
# Only run specific step-1 databases
python main.py --step 1 --databases biogrid,stitch,lincs,gem_dist1

# Limit to N metabolites (for quick test runs)
python main.py --step 1 --max-metabolites 500

# Resume after an interrupted API run (skips completed caches)
python main.py --step 1 --resume

# Override number of worker threads
python main.py --all --workers 8

# Write output to a custom directory
python main.py --all --output-dir PKN_test
```

Available database keys for `--databases`:

```
biogrid  stitch  uniprot  intact  chembl  lincs  gem_dist1  gem_dist2  metalinks
```

---

## Database Sources

### Step 1 – Metabolite–Gene Interactions

#### BioGRID (`biogrid.py`)
- **Input**: Local `BIOGRID-CHEMICALS-*.chemtab.txt`
- **Matching**: InChIKey (from HMDB) → BioGRID Chemical InChIKey
- **URL format**: `https://thebiogrid.org/chemical/{BioGRID_ID}/{Chemical_Name}.html`

#### STITCH (`stitch.py`)
- **Input**: Local `9606.protein_chemical.links.v5.0.tsv`
- **Matching**: PubChem CID (from HMDB) → STITCH chemical ID; Ensembl protein IDs mapped to gene symbols via MyGene.info API
- **URL format**: `http://stitch.embl.de/cgi/show_network_section.pl?identifier=CID{pubchem_id}&species=9606`

#### UniProtKB (`uniprot.py`)
- **Input**: UniProtKB REST API (`/uniprotkb?query=inchikey:{key}&organism_id=9606`)
- **Matching**: InChIKey per metabolite; returns reviewed human proteins
- **URL format**: `https://www.uniprot.org/uniprotkb/{accession}/entry` (per protein)

#### IntAct (`intact.py`)
- **Input**: EBI IntAct REST API (`/intact/ws/interaction/findInteractorsByInteractionIdentifier/CHEBI:{id}`)
- **Matching**: ChEBI ID per metabolite; extracts interacting human gene symbols and EBI interaction accession IDs
- **URL format**: `https://www.ebi.ac.uk/intact/details/interaction/{EBI_ID}` (per interaction)

#### ChEMBL (`chembl.py`)
- **Input**: ChEMBL REST API (bioactivity endpoint, filtered by organism = *Homo sapiens*)
- **Matching**: ChEMBL compound ID (mapped from canonical SMILES in preprocessing)
- **URL format**: `https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}`

#### LINCS (`lincs.py`)
- **Input**: Local `lsp_biochem_agg.csv`, `lsp_compound_mapping.csv`, `lsp_target_dictionary.csv`
- **Matching**: ChEMBL ID → LINCS compound ID (`lspci_id`) → gene targets below IC50 threshold (default 10 µM); human targets only
- **URL format**: `https://labsyspharm.github.io/lsp-cheminformatics/compounds/{lspci_id}`

#### Human-GEM (`gem.py`)
- **Input**: Local Human-GEM model files (`Human-GEM.txt`, `metabolites.tsv`, `genes.tsv`)
- **Method**: Builds a NetworkX directed graph of the metabolic network; for each metabolite finds all genes reachable within the configured pathway distance (1 or 2 reaction hops) via ego-graph traversal; ubiquitous currency metabolites (ATP, NAD+, etc.) are excluded from the graph
- **URL**: `https://doi.org/10.1126/scisignal.aaz1482` (paper) — same for all GEM edges
- **Extra column**: `reaction_id` — semicolon-separated list of reaction IDs connecting metabolite to gene

#### MetalinksDB (`metalinks.py`)
- **Input**: Local `metalinks.csv` (from the liana+ R package)
- **Matching**: HMDB ID directly; the `source` column in metalinks.csv is used to assign per-source provenance URLs
- **URL mapping**:
  - `CellPhoneDB` → `https://www.cellphonedb.org`
  - `NicheNet` → `https://github.com/saeyslab/nichenetr`
  - `HMDB` → `https://hmdb.ca`
  - others → `https://metalinks.omnipathdb.org`

### Step 2 – Protein–Protein Interactions

#### BioGRID PPI (`biogrid_ppi.py`)
- **Input**: Local `BIOGRID-ALL-*.tab3.txt`, filtered to *Homo sapiens* physical interactions
- **Output**: Gene symbol pairs with `Source = BioGRID_PPI`

#### HuRI (`huri.py`)
- **Input**: Local `HuRI_Tong_2021.tsv`; Ensembl IDs mapped to HGNC symbols via local `ensembl_to_hgnc.tsv`
- **Output**: Gene symbol pairs with `Source = HuRI`

#### HumanNet (`humannet.py`)
- **Input**: Local `HS-PI.symbol.tsv` (already in gene symbol format)
- **Output**: Gene symbol pairs with `Source = HumanNet`

#### STRING (`string_api.py`)
- **Input**: STRING v12 REST API; gene list derived from step-1 output
- **Filter**: Minimum combined confidence score 400 (medium confidence); human only (`species=9606`)
- **Output**: Gene symbol pairs with `Source = STRING`, `combined_score` column

### Step 3 – Assembly

#### Combiner (`combiner.py`)
Reads `metabolite_gene_PKN.tsv` (step 1) and `PPI_network.tsv` (step 2), renames columns to the standardized `Node1 / Node2 / Type / Source` format, concatenates them preserving optional columns (`url`, `reaction_id`, `combined_score`), and writes `LemonIte_PKN.tsv`.

#### Annotator (`annotator.py`)
Loads `LemonIte_PKN.tsv`, renames the lowercase `url` column to `URL`, fills missing values, logs per-source URL coverage statistics, and writes the final `LemonIte_PKN_with_URLs.tsv`. **No URL regeneration is performed here** — all URLs originate from the individual retrievers in step 1.

---

## Output Files

All files are written to `$PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME/` (default: `…/PKN/`).

### Intermediate (step 1)

| File | Columns | Description |
|---|---|---|
| `HMDB_metabolites_<DB>_processed.csv` | `HMDB_ID, Gene, Source, url [, reaction_id]` | Per-database cache; regenerated if absent |
| `HMDB_metabolites_<DB>_processed.csv.meta` | JSON | Cache metadata (timestamp, row count) |
| `metabolite_gene_PKN.tsv` | `Metabolite, Gene, Source, url, reaction_id` | Merged step-1 output |

### Intermediate (step 2)

| File | Columns | Description |
|---|---|---|
| `PPI_network.tsv` | `GeneA, GeneB, Source [, combined_score]` | Merged PPI network |
| `cache/string_ppi_cache.csv` | — | STRING API result cache |

### Final outputs (step 3)

| File | Columns | Description |
|---|---|---|
| `LemonIte_PKN.tsv` | `Node1, Node2, Type, Source, url, reaction_id, combined_score` | Combined PKN before annotation |
| `LemonIte_PKN_with_URLs.tsv` | `Node1, Node2, Type, Source, URL, reaction_id, combined_score` | **Primary output** — ready for modelling |

### Log

| File | Description |
|---|---|
| `pkn_pipeline.log` | Full pipeline log with timestamps |

---

## Caching System

Each step-1 retriever caches its result to `HMDB_metabolites_<DB>_processed.csv`. On the next run the cache is loaded directly, skipping the database query. A `.meta` JSON sidecar records the creation timestamp and row count.

**Cache invalidation** — delete the cache file (and its `.meta`) to force regeneration:

```bash
# Force a single database to re-query
rm PKN/HMDB_metabolites_UniProtKB_processed.csv
rm PKN/HMDB_metabolites_UniProtKB_processed.csv.meta

# Force all step-1 databases (keeps PPI and final outputs)
rm PKN/HMDB_metabolites_*_processed.csv*
```

> **Important**: caches generated before URL support was added (columns `HMDB_ID, Gene, Source` only) must be deleted. Current caches include `url` (and `reaction_id` for GEM) as required columns.

---

## HPC Deployment

### Standard run (full PKN)

```bash
cd PKN_build_pipeline
qsub submit_pkn_hpc.sh
```

The job requests 1 node × 24 cores, 128 GB RAM, 72 h walltime on the VSC Ghent (Doduo) cluster.

Key environment variables set by `submit_pkn_hpc.sh`:

```bash
PKN_WORKDIR=/user/gent/435/vsc43501/25BVM_Lemonite
PKN_OUTPUT_DIR_NAME=PKN
PKN_DB_DIR=$LEMONITE_DIR/databases
PKN_GEM_DIR=$LEMONITE_DIR/models/Human1-GEM/model
PKN_CONFIG=config_hpc
```

### Test run (subset of metabolites)

A separate submission script runs step 1 on the first 2000 metabolites only:

```bash
qsub submit_pkn_hpc_2000.sh   # output → PKN_test_2000/
```

### Monitoring

```bash
qstat -u $USER              # check job status
tail -f PKN/pkn_pipeline.log  # live log
```

---

## Singularity Container

A `Singularity.def` recipe is included for fully portable, reproducible execution.

### Build

```bash
singularity build pkn_pipeline.sif Singularity.def
```

### Run

```bash
singularity run \
    --bind /path/to/databases:/databases \
    --bind /path/to/output:/output \
    pkn_pipeline.sif --all
```

Environment variables can be passed with `--env`:

```bash
singularity run \
    --bind /path/to/databases:/databases \
    --bind /path/to/output:/output \
    --env PKN_WORKDIR=/output \
    --env PKN_DB_DIR=/databases \
    pkn_pipeline.sif --all
```

---

## Troubleshooting

### Empty URL column in output

The most common cause is a stale cache missing the `url` column. Delete the relevant cache files:

```bash
rm PKN/HMDB_metabolites_<DB>_processed.csv*
```

Then re-run the pipeline; the retriever will regenerate the cache with proper URLs.

### API timeouts / connection errors

The `@retry_api_call` decorator retries automatically with exponential backoff. If an API consistently fails, reduce `max_workers` or increase `backoff_factor` in `config.py`. Use `--resume` to skip databases whose caches are already complete.

### HMDB XML not found

Ensure `HMDB_METABOLITES_XML` (or `PKN_DB_DIR`) points to the correct location. The HMDB XML can be downloaded from [hmdb.ca/downloads](https://hmdb.ca/downloads).

### GEM graph build fails

Check that `PKN_GEM_DIR` (or `GEM_PATH`, `GEM_METABOLITES_PATH`, `GEM_GENES_PATH` in config) points to the Human-GEM v1.x model directory containing `Human-GEM.txt`, `metabolites.tsv`, and `genes.tsv`.

### Step 2 fails ("metabolite-gene PKN not found")

Step 2 requires step 1 output. Run `python main.py --step 1` first, or use `python main.py --all`.

### Step 3 fails ("PPI network not found")

Step 3 requires both step-1 and step-2 output. Run steps in order, or use `python main.py --all`.

---

## Repository Layout

```
build_PKN/
├── PKN_build_pipeline/       ← This pipeline
│   ├── README.md             ← This file
│   ├── QUICKSTART.md         ← Brief getting-started guide
│   ├── main.py
│   ├── config.py
│   ├── config_hpc.py
│   ├── requirements.txt
│   ├── setup.sh
│   ├── Singularity.def
│   ├── submit_pkn_hpc.sh     ← Full HPC job (all metabolites)
│   ├── submit_pkn_hpc_2000.sh← Test HPC job (first 2000 metabolites)
│   ├── step1_metabolites/
│   ├── step2_proteins/
│   ├── step3_final/
│   └── utils/
│
├── Collect_PKNdata_metabolites.ipynb  ← Reference notebook (step 1)
├── Collect_PKNdata_proteins.ipynb     ← Reference notebook (step 2)
├── Build_final_PKN.ipynb              ← Reference notebook (step 3)
└── test_minimal.py
```

The Jupyter notebooks in the parent directory were the original exploratory implementations and serve as the reference for URL formats and retrieval logic. The pipeline code in `PKN_build_pipeline/` is the production version that produces identical results.
