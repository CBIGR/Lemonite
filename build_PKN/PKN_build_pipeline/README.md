# PKN Pipeline - Modularized Python Scripts

## Overview

This directory contains a modularized Python implementation of the 3-step PKN (Prior Knowledge Network) pipeline, refactored from the original Jupyter notebooks.

## Current Implementation Status

### ✅ Phase 1: Infrastructure (COMPLETED)
- [x] Directory structure created
- [x] `config.py` - Centralized configuration
- [x] `utils/api_retry.py` - API retry decorators
- [x] `utils/pipeline.py` - Base retriever classes
- [x] `utils/file_io.py` - File I/O utilities
- [x] `main.py` - CLI orchestrator
- [x] `requirements.txt` - Python dependencies

### 🚧 Phase 2: Step 1 Retrievers (IN PROGRESS)
- [x] `step1_metabolites/biogrid.py` - BioGRID retriever
- [ ] `step1_metabolites/stitch.py` - STITCH retriever
- [ ] `step1_metabolites/uniprot.py` - UniProt retriever
- [ ] `step1_metabolites/intact.py` - IntAct retriever
- [ ] `step1_metabolites/chembl.py` - chEMBL retriever
- [ ] `step1_metabolites/lincs.py` - LINCS retriever
- [ ] `step1_metabolites/l1000.py` - L1000 gene expression retriever
- [ ] `step1_metabolites/gem.py` - Human-GEM pathway retriever
- [ ] `step1_metabolites/metalinks.py` - MetalinksDB retriever
- [ ] `step1_metabolites/preprocessing.py` - Metabolite preprocessing
- [ ] `step1_metabolites/integration.py` - Database integration

### ⏳ Phase 3: Step 2 PPI Retrievers (PENDING)
- [ ] `step2_proteins/string_api.py` - STRING API retriever
- [ ] `step2_proteins/biogrid_ppi.py` - BioGRID PPI retriever
- [ ] `step2_proteins/huri.py` - HuRI retriever
- [ ] `step2_proteins/ppi_integration.py` - PPI integration

### ⏳ Phase 4: Step 3 Final PKN (PENDING)
- [ ] `step3_final/combiner.py` - Network combiner
- [ ] `step3_final/annotator.py` - URL annotator
- [ ] `step3_final/analysis.py` - Comparative analysis
- [ ] `step3_final/visualization.py` - Visualization

## Installation

```bash
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Run all steps
```bash
python main.py --all
```

### Run individual steps
```bash
# Step 1: Metabolite-gene interactions
python main.py --step 1

# Step 2: Protein-protein interactions
python main.py --step 2

# Step 3: Final PKN integration
python main.py --step 3
```

### Run specific databases only
```bash
python main.py --step 1 --databases biogrid,stitch,lincs
```

### Resume from checkpoint (for L1000)
```bash
python main.py --step 1 --resume
```

## Architecture

### Infrastructure Layer
```
config.py                    # All paths and settings
utils/
  ├── api_retry.py          # Retry decorators with exponential backoff
  ├── pipeline.py           # Base classes (DatabaseRetriever, LocalFileRetriever, APIRetriever)
  └── file_io.py            # HMDB loading, result saving, progress tracking
```

### Step 1: Metabolite-Gene Interactions
```
step1_metabolites/
  ├── biogrid.py            # BioGRID CHEMICALS
  ├── stitch.py             # STITCH protein-chemical links
  ├── uniprot.py            # UniProt XML + API
  ├── intact.py             # IntAct API
  ├── chembl.py             # chEMBL API
  ├── lincs.py              # LINCS biochemical data
  ├── l1000.py              # L1000 gene expression (2.1GB GMT)
  ├── gem.py                # Human-GEM pathway distances
  ├── metalinks.py          # MetalinksDB
  ├── preprocessing.py      # ChEMBL ID mapping, SMILES canonicalization
  └── integration.py        # Combine all databases into final PKN
```

### Step 2: Protein-Protein Interactions
```
step2_proteins/
  ├── string_api.py         # STRING API (chunks of 1000 genes)
  ├── biogrid_ppi.py        # BioGRID PPI local file
  ├── huri.py               # HuRI with Ensembl mapping
  └── ppi_integration.py    # Combine PPI sources + Venn diagram
```

### Step 3: Final PKN Integration
```
step3_final/
  ├── combiner.py           # Vertical concatenation of metabolite-gene + PPI
  ├── annotator.py          # Add IntAct/STITCH/BioGRID/UniProt/chEMBL URLs
  ├── analysis.py           # Compare vs MetalinksDB/MEBOCOST
  └── visualization.py      # Superclass heatmaps, network stats
```

## Configuration

Edit `config.py` to change:
- Output directory: `OUTPUT_DIR_NAME = 'PKN'`
- Database file paths
- API retry settings
- Thread pool sizes
- Logging levels

## Key Features

### Retry Logic
All API retrievers use exponential backoff with configurable:
- Max retries
- Backoff factor
- Timeout duration
- Rate limiting

### Caching
- Local file retrievers cache parsed results
- API retrievers cache raw responses
- Resume capability for long-running operations (L1000)

### Base Classes
All retrievers inherit from:
- `LocalFileRetriever`: For local files (BioGRID, STITCH, GEM)
- `APIRetriever`: For REST APIs (UniProt, IntAct, chEMBL, STRING)

Both implement:
- `get_interactions(metabolites) -> pd.DataFrame`
- Automatic caching via `load_cache()` / `save_cache()`

## Output Files

### Step 1 Outputs
```
PKN/
  ├── HMDB_metabolites_BioGRID_processed.csv
  ├── HMDB_metabolites_STITCH_processed.csv
  ├── HMDB_metabolites_UniProtKB_processed.csv
  ├── HMDB_metabolites_IntAct_processed.csv
  ├── HMDB_metabolites_chEMBL_processed.csv
  ├── HMDB_metabolites_LINCS_processed.csv
  ├── HMDB_metabolites_L1000_processed.csv
  ├── HMDB_metabolites_Human1_GEM_dist1_processed.csv
  ├── HMDB_metabolites_Human1_GEM_dist2_processed.csv
  ├── HMDB_metabolites_MetalinksDB_processed.csv
  └── metabolite_gene_PKN.tsv
```

### Step 2 Outputs
```
PKN/
  └── PPI_network.tsv
```

### Step 3 Outputs
```
PKN/
  ├── LemonIte_PKN.tsv
  └── LemonIte_PKN_with_URLs.tsv
```

## Logging

Logs are written to:
- `PKN/api_errors.log` - API failures and retries
- `PKN/pipeline_progress.log` - High-level pipeline progress

## Development Workflow

### Adding a New Database Retriever

1. Create new file in `step1_metabolites/` (e.g., `newdb.py`)
2. Inherit from `LocalFileRetriever` or `APIRetriever`
3. Implement `get_interactions(metabolites) -> pd.DataFrame`
4. Add to `main.py` in the `run_step1_metabolites()` function
5. Test standalone: `python step1_metabolites/newdb.py`

### Testing Individual Retrievers

Each retriever has a `if __name__ == '__main__':` block for standalone testing:

```bash
python step1_metabolites/biogrid.py
```

## Notes

- **L1000**: 2.1GB GMT file, single-threaded to avoid RAM issues
- **STITCH**: Requires MyGene.info for Ensembl→gene symbol mapping
- **GEM**: NetworkX graph with pathway distance calculations
- **APIs**: UniProt, IntAct, chEMBL, STRING - all have retry logic

## Next Steps

1. Complete Step 1 retrievers (9 remaining)
2. Implement preprocessing and integration modules
3. Build Step 2 PPI retrievers
4. Build Step 3 final PKN modules
5. Comprehensive testing
6. Documentation and examples

## License

See parent project for license information.
