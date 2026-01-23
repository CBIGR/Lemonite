# Lemonite Pipeline - Complete Documentation

> **Note**: This is the full documentation for the Lemonite pipeline. For quick start instructions, see the [main README](README.md).

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Installation & Setup](#installation--setup)
3. [Input Data Formats](#input-data-formats)
4. [Output Structure](#output-structure)
5. [Pipeline Parameters](#pipeline-parameters)
6. [Advanced Usage](#advanced-usage)
7. [Troubleshooting](#troubleshooting)
8. [Performance Optimization](#performance-optimization)
9. [Container Management](#container-management)

---

## Pipeline Overview

### What is Lemonite?

Lemonite is a comprehensive Nextflow pipeline for multi-omics data integration that identifies gene co-expression modules and their mulit-omics regulators. Building upon the [LemonTree algorithm](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003983) (Bonnet et al., 2015), Lemonite extends it with:

- **Multi-regulator support**: Dynamically handle TFs, metabolites, lipids, and custom regulator types
- **Flexible regulator selection**: Multiple strategies for regulator selection (percentage-based, fold-change per module)
- **Comprehensive evaluation**: Validation against the Lemonite prior knowledge networks (PKNs), specifically for metabolomics data
- **Functional enrichment analyses**: Pathway analysis on module genes and regulator target genes
- **Rich visualizations**: Module heatmaps, network graphs, interactive network figures

### Detailed Pipeline Workflow

```
📊 Input Data (Transcriptomics + Metabolomics + Optional Omics)
    ↓
1. PREPROCESSING & TFA
    • DESeq2 normalization (transcriptomics)
    • Log transformation with auto-detection for already-logged data
    • Omics-specific scaling (Pareto for metabolomics/lipids, Z-score for transcriptomics)
    • Transcription factor activity inference using decoupleR (optional, with CollecTRI or custom network)
    • Highly variable gene selection (top N genes by variance)
    • Vertical integration of all omics datasets into single matrix
    ↓
2. LEMONTREE CLUSTERING
    • Parallel Gibbs sampler for simultaneous inference of condition and co-expression clusters (10-100 independent runs)
    • Consensus clustering across all runs
    • Consensus regulator assignment for each regulator type (TFs, metabolites, etc.) using an ensembl of decision trees
    • Generation of real and random regulator scores for each co-expression module
    ↓
3. NETWORK GENERATION & FILTERING
    • Filtering out of low-quality modules (default coherence threshold: 0.6)
    • Regulator selection using configurable method:
      - Percentage-based: Top N% of regulators
      - Fold-change per module: Score ≥ (max_random(module) × fold_cutoff) per module
    • Z-score normalization of selected regulator scores:
      - Makes scores more comparable across different omics types
      - Uses (score - mean) / std transformation
    • Network construction with regulator→gene edges
    • Hub analysis to identify highly connected regulators
    • Export to Cytoscape-compatible formats and files used in downstream analyses
    ↓
4. DOWNSTREAM ANALYSES
    ├─→ PKN EVALUATION
    │   • Validate metabolite-gene connections against Lemonite knowledge graph
    │   • Generate files for module visualizations
    │
    ├─→ SUBNETWORK GRAPHS
    │   • Visualize predicted regulator-module gene connections in the Lemonite knowledge graph
    │   • Edge categorization (Causal, Metabolic pathway, PPI, Ambiguous)
    │   • Color-coded nodes by regulator type
    │   • Per-module Cytoscape files
    │
    ├─→ MODULE HEATMAPS
    │   • Gene expression and regulator expression heatmaps per module
    │   • Sample ordering by eigengene
    │   • Annotation tracks (sample groups, regulators)
    │   • Optional metabolite-gene interaction overlay
    │
    └─→ ENRICHMENT ANALYSIS
        • Pathway enrichment on module genes (up/down regulated)
        • Regulator target gene set enrichment
        • Multiple databases: GO, KEGG, Reactome
        • EnrichR and/or GSEA methods available
    ↓
5. MODULE OVERVIEW
    • Interactive HTML visualization of complete integrated regulatory network
    • Module prioritization based on differential expression (if possible)
    • Functional similarity-based module clustering (optional megago)
    • Summary table with all module statistics
    • Expression heatmap across all modules
```

---

## Installation & Setup

### System Requirements

**Minimum:**
- 16 GB RAM
- 4 CPU cores

**Recommended:**
- 32 GB RAM
- 8+ CPU cores

### Software Prerequisites

1. **Nextflow** (≥ 25.04.0)
2. **Container runtime**: Docker OR Singularity
3. **Java** 11+ (for Nextflow)

### Step-by-Step Installation

#### 1. Install Nextflow

```bash
# Download and install Nextflow
curl -s https://get.nextflow.io | bash

# Make executable and move to PATH
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

#### 2. Install Container Runtime

**Docker (Local/Workstation):**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install docker.io
sudo systemctl start docker
sudo systemctl enable docker

# Add user to docker group (avoid sudo)
sudo usermod -aG docker $USER
# Log out and back in for group changes to take effect

# Verify
docker --version
```

**Singularity (HPC):**
```bash
# Usually pre-installed on HPC clusters
# Check with your system administrator
singularity --version

# If not available, see: https://docs.sylabs.io/guides/latest/admin-guide/
```

#### 3. Clone Repository

```bash
git clone https://github.com/CBIGR/Lemonite.git
cd Lemonite/nextflow
```

#### 4. Build Container Images

**Docker:**
```bash
./build-docker.sh
# Or manually:
docker build -t lemontree-pipeline:latest .
```

**Singularity:**
```bash
./build-singularity.sh
# Or manually:
singularity build lemontree-pipeline.sif Singularity.def
```

---

## Input Data Formats

### Directory Structure

Your input directory must follow this structure:

```
input_directory/
└── data/
    ├── Counts.tsv               # REQUIRED (data for clustering into co-expression modules)
    ├── Metabolomics.txt         # OPTIONAL (if using metabolites as regulators)
    ├── Lipidomics.txt           # OPTIONAL (if using lipids as regulators)
    ├── Proteomics.txt           # OPTIONAL (if using proteins as regulators)
    ├── [Other omics files]      # OPTIONAL (any additional omics types)
    ├── Metadata.txt             # REQUIRED (sample metadata)
    ├── Lovering_TF_list.txt     # OPTIONAL (List if TFs to assign as regulator, default is list from [Lovering_TF_list.txt](https://www.sciencedirect.com/science/article/pii/S1874939921000833?via%3Dihub))
    └── name_map.csv             # OPTIONAL (Required for comparison of metabolite regulators with Lemonite knowledge graph)
```

### Detailed File Specifications

#### 1. Counts.tsv (Gene Expression - protein abundance)

**Format:** Tab-separated values  
**Rows:** Genes (HGNC symbols)  
**Columns:** Sample IDs  
**Values:** Raw read counts (integers) or protein abundances

**Important:**
- First column contains gene symbols (no column name)
- Column headers must match Metadata sample IDs (and also other omics if these are provided)
- Ideally raw counts (not TPM, RPKM, or FPKM)
- Remove genes with all zeros

**Example:**
```tsv
	Sample1	Sample2	Sample3	Sample4
TP53	1234	2345	3456	4567
EGFR	567	678	789	890
MYC	2345	3456	4567	5678
STAT3	890	901	1012	1123
```


#### 2. Metabolomics.txt (Metabolite Abundance)

**Format:** Tab-separated values  
**Rows:** Metabolites (feature names will be cleaned by replacing '(),-:\/' by '_')
**Columns:** Sample IDs (matching Counts.tsv and Metadata.txt)  
**Values:** Abundances (raw intensities or log-transformed, NAs will be replaced by 0s)

**Important:**
- First column can have header "Name" or will be treated as metabolite names
- Sample IDs must exactly match gene expression and metadata files
- Can be raw intensities OR already log-transformed
- Pipeline auto-detects log-transformation (median < 100 → already logged)

**Example (raw intensities):**
```tsv
Name	Sample1	Sample2	Sample3	Sample4
L-glutamine	1234567	2345678	3456789	4567890
Citric_acid	567890	678901	789012	890123
Pyruvate	234567	345678	456789	567890
```

**Example (log2-transformed):**
```tsv
Name	Sample1	Sample2	Sample3	Sample4
L-glutamine	20.23	21.16	21.72	22.13
Citric_acid	19.12	19.37	19.60	19.77
Pyruvate	17.84	18.40	18.80	19.12
```

**Auto-detection logic:**
```python
if median(data) < 100 and max(data) < 100:
    # Data is already log-transformed, skip log transformation
else:
    # Apply log2(x + 1) transformation
```

#### 3. Additional Omics Files (Optional)

**Format:** Same as Metabolomics.txt (tab-separated, features as rows, samples as columns) - samples IDs must match 
**Purpose:** Additional omics to assign as regulators
**Note:** Must be continuous data (not binary)
**File naming:** Can use any descriptive name (e.g., `Proteomics.txt`, `Lipidomics.tsv`, `Kinases.txt`)
**Configuration:** Specify in `--regulator_types` parameter to include in analysis

#### 4. Metadata.txt (Sample Metadata)

**Format:** Tab-separated values  
**Rows:** Samples (matching expression/metabolomics)  
**Columns:** Sample attributes

**Required columns:**
- `Sample_ID`: Sample identifiers (matching data files)
- `diagnosis` (or custom condition column): Primary grouping variable

**Optional columns:**
- `age`, `sex`, `batch`, `biopsy_location`, etc.
- Used for stratification, visualization, and DESeq2 analysis
- Column names for metadata columns to include can be specified via `--metadata_columns` parameter
- These columns will be saved in `DESeq_groups.txt` and `sample_mapping.mvf` for downstream visualization

**Important Note:**
- The `--metadata_columns` parameter controls which columns are **saved** from your metadata file during preprocessing
- To **display** specific metadata columns as annotation bars in module heatmaps, use the module_viewer.py `--annotation_types` parameter (see [Module Heatmap Annotations](#module-heatmap-annotations))
- These are independent: you can save many metadata columns but choose to display only some in visualizations

**Example:**
```tsv
Sample_ID	diagnosis	age	sex	batch	biopsy_location
Sample1	case	45	M	1	colon
Sample2	control	38	F	1	colon
Sample3	case	52	M	2	rectum
Sample4	control	41	F	2	rectum
```

**DESeq2 Design Formula:**
- Default: `~ diagnosis`
- Custom: `--design_formula "~ batch + diagnosis"`
- Complex: `--design_formula "~ biopsy_location + diagnosis + sex"`

#### 5. lovering_TF_list.txt (Transcription Factor List)

**Format:** Plain text, one TF per line, can be a custom list or default Lovering list  
**Content:** HGNC gene symbols of transcription factors

**Example:**
```txt
TP53
MYC
STAT3
NFKB1
JUN
FOS
```

**Custom TF lists:**
```bash
# Use your own list
cp my_custom_TFs.txt data/MyTFs_list.txt
--regulator_types "TFs:MyTFs_list.txt,Metabolites:Metabolomics.txt"
```

**Default list:** [Lovering_TF_list.txt](https://www.sciencedirect.com/science/article/pii/S1874939921000833?via%3Dihub) - curated human TF list  
**Location:** `PKN/lovering_TF_list.txt` (1,456 TFs)

#### 7. name_map.csv (Metabolite ID Mapping)

**Format:** CSV  
**Purpose:** Map metabolite names to HMDB IDs for PKN validation

**Required columns:**
- `Query`: Metabolite name (as in Metabolomics.txt)
- `HMDB`: HMDB identifier

**Example:**
```csv
Query,HMDB
L-glutamine,HMDB0000641
Citric_acid,HMDB0000094
Pyruvate,HMDB0000243
Glucose,HMDB0000122
```

**How to generate:**
1. Go to [MetaboAnalyst Compound ID Conversion](https://www.metaboanalyst.ca/MetaboAnalyst/upload/ConvertView.xhtml)
2. Upload metabolite list (from Metabolomics.txt first column)
3. Select "Name" as input type
4. Generate name mapping, make sure to mannualy refine after
5. Download CSV result and save as name_map.csv in data directory


---

## Pipeline Parameters

### Complete Parameter Reference

#### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--input_dir` | Path | **REQUIRED** | Directory containing `data/` folder with input files |
| `--output_dir` | Path | `./results` | Output directory for all results |
| `--data_dir` | Path | `{input_dir}/data` | Specific data directory (auto-detected if not provided) |
| `--run_id` | String | Auto-generated | Unique identifier for this run (used in output naming) |
| `--lemontree_jar` | Path | `/opt/lemontree/lemontree_v3.1.1.jar` | Path to LemonTree JAR file |
| `--max_cpus` | Integer | `16` | Maximum CPUs for parallel tasks |
| `--max_memory` | String | `64.GB` | Maximum memory allocation |
| `--max_time` | String | `24.h` | Maximum time for any single task |
| `--organism` | String | `human` | Currently either human or mouse |

Run_id will be auto-generated if not specified by the user:
```
Format: {top_n_genes}HVG_coherence{threshold}_{method}_{n_clusters}
Example: 1500HVG_coherence0.6_fold2.5x_permodule_100
```

### Organism parameter

- **Effect on preprocessing**: selects the appropriate Ensembl dataset (`hsapiens_gene_ensembl` or `mmusculus_gene_ensembl`) and gene symbol attribute mapping when converting gene symbols to Ensembl IDs. If TFA is enabled with default prior network, will use the mouse version of CollecTri.
- **Effect on enrichment**: selects the appropriate organism-specific OrgDb (e.g., `org.Hs.eg.db` or `org.Mm.eg.db`) and uses the correct KEGG organism code for GSEA (`hsa` or `mmu`). EnrichR database selection will prefer GO databases for non-human species.
- **Effect on TF list**: automatically selects the appropriate transcription factor list (`Lovering_TF_list.txt` for human, `lovering_TF_list_mouse.txt` for mouse). This ensures TF symbols match your expression data.

**Accepted values:**
- `human` (default): Uses human gene annotations, KEGG/Reactome pathways, and human TF list
- `mouse`: Uses mouse gene annotations, KEGG/Reactome pathways, and mouse TF list

Usage example:
```bash
# Human analysis (default)
nextflow run main.nf --input_dir /path/to/data --organism human

# Mouse analysis
nextflow run main.nf --input_dir /path/to/data --organism mouse
```

**Important for Mouse Data:**
- The pipeline automatically uses `lovering_TF_list_mouse.txt` (1,141 mouse TF symbols) when `--organism mouse` is specified
- This file is generated from human-to-mouse TF mapping and uses proper mouse gene symbols (e.g., `Adnp`, `Ahr`, `Alx1`)
- If you have mouse expression data but don't specify `--organism mouse`, the pipeline will use human TF symbols and find no matches, causing errors
- The mouse TF list file must be present in your `data/` directory (automatically copied from `PKN/lovering_TF_list_mouse.txt` if missing)

#### Regulator Configuration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--regulator_types` | String | `"TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"` | Comma-separated Prefix:DataFile pairs |

**Format Specification:**
```
Prefix1:DataFile1,Prefix2:DataFile2,Prefix3:DataFile3,...
```

Where:
- **Prefix**: Used throughout pipeline for naming intermediate/output files (e.g., `TFs`, `Metabolites`, `Lipids`, `Proteins`)
- **DataFile**: Name of the input data file in your `data/` directory

**Important:**
- For **TFs**: DataFile should be a list file (e.g., `Lovering_TF_list.txt`)
- For **other regulators**: DataFile should be abundance data (e.g., `Metabolomics.txt`, `Proteomics.txt`)

**Examples:**

```bash
# Default: TFs and metabolites
--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"

# Add lipidomics
--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Lipids:Lipidomics.txt"

# Add proteomics
--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Proteins:Proteomics.txt"

# All omics types
--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Lipids:Lipidomics.txt,Proteins:Proteomics.txt"

# Custom TF list
--regulator_types "TFs:MyCustomTFs.txt,Metabolites:Metabolomics.txt"

# Single regulator type (TFs only)
--regulator_types "TFs:Lovering_TF_list.txt"
```

**File mapping example:**

For `Metabolites:Metabolomics.txt`:
- Input file: `data/Metabolomics.txt` (abundance data)
- Preprocessing creates: `Preprocessing/metabolites.txt` (regulator list extracted from row names)
- LemonTree creates: `Lemon_out/Metabolites.allreg.txt`, `Lemon_out/Metabolites.randomreg.txt`
- Network files: `Metabolites2targets_*.txt`

For `TFs:Lovering_TF_list.txt`:
- Input file: `data/Lovering_TF_list.txt` (list of TF gene symbols)
- Preprocessing copies to: `Preprocessing/lovering_TFs.txt`
- LemonTree creates: `Lemon_out/TFs.allreg.txt`, `Lemon_out/TFs.randomreg.txt`
- Network files: `TFs2targets_*.txt`

#### LemonTree Clustering Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--n_clusters` | Integer | `100` | Number of parallel Gibbs sampling runs |
| `--coherence_threshold` | Float | `0.6` | Minimum module coherence score (0-1) |

**n_clusters guidance:**
- Minimum: 10 runs (fast, less robust, test purposes only)
- Recommended: 100 runs (balance of speed/quality)
- High quality: 200+ runs (slow, best consensus)

**coherence_threshold:**
- Range: 0.0 (all modules) to 1.0 (only perfect co-expression)
- Recommended: 0.4-0.7
- Lower = more modules, potentially noisier
- Higher = fewer modules, higher quality

#### Preprocessing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--top_n_genes` | Integer | `5000` | Number of highly variable genes to select |
| `--use_omics_specific_scaling` | Boolean | `true` | Use Pareto for metabolomics and other omics, Z-score for transcriptomics |
| `--perform_tfa` | Boolean | `true` | Perform TF activity inference |
| `--deseq_contrast1` | String | `diagnosis` | Main metadata column for DESeq2 |
| `--design_formula` | String | `~ diagnosis` | DESeq2 design formula |
| `--metadata_columns` | String | `diagnosis` | Comma-separated metadata columns to include |
| `--expression_col` | String | `count` | Gene symbol column name in expression file |
| `--sample_id_col` | String | `Sample_ID` | Sample ID column in metadata |

**Omics-Specific Scaling Details:**

When `--use_omics_specific_scaling true` (recommended):

| Data Type | Normalization | Scaling | Formula | Rationale |
|-----------|---------------|---------|---------|-----------|
| Transcriptomics | DESeq2 → log2 | Z-score | (x - μ) / σ | Standard for RNA-seq |
| Metabolomics | log2(x+1) | Pareto | (x - μ) / √σ | Preserves variance structure |
| Lipidomics | log2(x+1) | Pareto | (x - μ) / √σ | Similar to metabolomics |
| Other | log2(x+1) | Pareto | (x - μ) / √σ | Default for unknowns |

**Why Pareto scaling for metabolomics?**
1. Less aggressive than z-score, preserves biological variance
2. Reduces influence of outliers (common in metabolomics)
3. Maintains relative magnitude differences between metabolites
4. Standard practice in metabolomics field
5. Better for data with heteroscedastic variance

When `--use_omics_specific_scaling false` (legacy):
- All omics types use Z-score scaling
- Not recommended for metabolomics/lipidomics data

**TFA (Transcription Factor Activity) Inference:**
- Uses decoupleR with consensus
- Default network: CollecTRI (curated TF-target interactions)
- Generates `TFA_consensus.txt` with inferred TF activities
- TFA values replace gene expression for TFs in clustering

**DESeq2 Design Formula Examples:**
```bash
# Simple (default)
--design_formula "~ diagnosis"

# With batch correction
--design_formula "~ batch + diagnosis"

# Complex design
--design_formula "~ biopsy_location + diagnosis + sex"

# Interaction term
--design_formula "~ diagnosis * sex"
```

#### Regulator Selection Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--regulator_selection_method` | String | `percentage` | Selection strategy: `percentage` or `fold_per_module` |
| `--top_n_percent_regulators` | Float | `2.0` | Top N% for percentage method |
| `--regulator_fold_cutoff` | Float | `2.0` | Fold cutoff for fold_per_module method |

**Method Comparison:**

##### 1. Percentage Method
```bash
--regulator_selection_method percentage
--top_n_percent_regulators 2.0
```

**Algorithm:**
```python
for each module:
    scores = all_regulator_scores_for_module
    threshold = percentile(scores, 100 - top_n_percent)
    selected = regulators where score >= threshold
```


##### 2. Fold-Change Per Module Method (Recommended)
```bash
--regulator_selection_method fold_per_module
--regulator_fold_cutoff 2.5
```

**Algorithm:**
```python
for each module:
    max_random = max(random_regulator_scores_for_module)
    threshold = max_random * regulator_fold_cutoff
    selected = regulators where score >= threshold
```

**Recommended cutoffs:**
- Conservative: 3.0-4.0 (fewer, high-confidence regulators)
- Moderate: 2.0-2.5 (balanced)
- Permissive: 1.5-2.0 (more regulators, some noise)

#### Z-Score Normalization of Regulator Scores

**Important:** Z-score normalization is applied **AFTER** regulator selection to preserve the original score distributions used for filtering.

**Purpose:**
- Makes regulator scores more comparable across different omics types (TFs, metabolites, lipids, proteins)
- Each omics type has different score ranges from LemonTree

**Algorithm:**
```python
# For each regulator type independently:
1. Select regulators using original LemonTree scores
2. Calculate mean and std from TRUE regulator scores only (exclude random scores)
3. Apply transformation: z_score = (score - mean) / std
4. Both true and random scores normalized using same parameters
```

**Z-score Interpretation:**
- Z = 0: Average regulator score
- Z = 1: One standard deviation above average (~84th percentile)
- Z = 2: Two standard deviations above average (~97.7th percentile)
- Z = 3: Three standard deviations above average (~99.9th percentile)

**Why Normalize After Selection:**
- Original scores reflect biological signals in each omics type
- Selection thresholds (percentile, fold-change) designed and tested for raw score distributions
- Post-selection normalization only affects visualization and comparability


**Files created:**
- `Lemon_out/{prefix}.selected_top{X}pct.txt`: Selected regulator-module pairs with normalized z-scores
- `ModuleViewer_files/{prefix}.selected_regs_list.txt`: List of selected regulator names
- Network edge files: `{prefix}2targets_{method}_{n_modules}_modules.txt`
- All downstream analyses use the normalized scores from the selected files

#### Enrichment Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--enrichment_method` | String | `EnrichR` | Method: `EnrichR`, `GSEA`, or `both` |

**EnrichR Databases (default):**
- GO_Biological_Process_2025
- GO_Molecular_Function_2025
- GO_Cellular_Component_2025
- KEGG_2021_Human
- Reactome_Pathways_2024

#### Module Overview Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--use_megago` | Boolean | `false` | Use [megago](http://megago.ugent.be/about) for functional clustering ⚠️ |
| `--prioritize_by_expression` | Boolean | `true` | Rank modules by differential expression |

**⚠️ Megago functional clustering:**
- **WARNING**: Computationally intensive! Can take 12-48+ hours for large datasets
- Only enable if you need advanced functional similarity analysis
- Groups modules by functional similarity (top 10 upregulated pathways)
- Uses GO term enrichment overlap while considering hierarchical organization of GO terms
- **For HPC systems**: Use `--use_megago true` with the `hpc` profile (see below)

#### PKN Evaluation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--pkn_network` | Path | `PKN/Lemonite_PKN.tsv` | Lemonite prior knowledge network |

---

## Output Structure

### Complete Output Directory Layout

```
results/
└── {run_id}/
    └── LemonTree/
        ├── Preprocessing/
        │   ├── LemonPreprocessed_expression_TFA.txt
        │   ├── LemonPreprocessed_metabolomics.txt
        │   ├── LemonPreprocessed_complete.txt
        │   ├── TFA_consensus.txt
        │   ├── Normalized_expression.pdf
        │   ├── Normalized_metabolomics.pdf
        │   ├── expression_variance_histogram.pdf
        │   ├── DESeq_groups.txt
        │
        ├── Lemon_out/
        │   ├── Lemon_results/
        │   │   ├── cluster_1/
        │   │   ├── cluster_2/
        │   │   └── ... (one per n_clusters)
        │   ├── Lovering.allreg.txt
        │   ├── Lovering.randomreg.txt
        │   ├── Metabolites.allreg.txt
        │   ├── Metabolites.randomreg.txt
        │   ├── tight_clusters.txt (consensus clusters)
        │   └── clusterfile
        │
        ├── Networks/
        │   ├── LemonNetwork_{method}_{N}modules.txt
        │   ├── TFs2targets_{method}_{N}_modules.txt
        │   ├── Metabolites2targets_{method}_{N}_modules.txt
        │   ├── Other_omics2targets_{method}_{N}_modules.txt
        │   ├── Cytoscape_edges_file.txt
        │   ├── Cytoscape_nodes_attributes_file.txt
        │   ├── specific_modules.txt
        │   ├── Module_coherence_scores.txt
        │   ├── *.selected_regs_*.png (most connected regulators per omics type)
        │   │
        │   └── subnetworks/
        │       ├── graph_module1_categorized_clean.png (visualization of module 1 in Lemonite knowledge graph)
        │       ├── graph_module2_categorized_clean.png
        │       └── cytoscape/
        │           ├── module_1_edges_categorized.tsv
        │           ├── module_1_attributes.tsv
        │           └── ...
        │
        ├── ModuleViewer_files/
        │   ├── clusters_list.txt
        │   ├── TFs.selected_regs_list.txt
        │   ├── Metabolites.selected_regs_list.txt
        │   ├── Lipids.selected_regs_list.txt
        │   ├── Proteins.selected_regs_list.txt
        │   ├── sample_mapping.mvf
        │   ├── evaluation_summary.txt
        │   └── metabolite_LemoniteKG_interactions.mvf
        │
        ├── module_heatmaps/
        │   └── heatmaps/
        │       ├── module_1_heatmap.png
        │       ├── module_2_heatmap.png
        │       └── ...
        │
        ├── Enrichment/
        │   ├── Modules_enrichr/
        │   │   ├── Enrichr_top_10_enriched_pathways_up_per_module.csv
        │   │   ├── module_1/
        │   │   │   ├── GO_Biological_Process_2025.png
        │   │   │   ├── GO_Molecular_Function_2025.png
        │   │   │   ├── GO_Cellular_Component_2025.png
        │   │   │   ├── KEGG_2021_Human.png
        │   │   │   └── Reactome_Pathways_2024.png
        │   │   └── ...
        │   ├── Modules_gsea/ (if enabled)
        │   │   ├── Enrichr_top_10_enriched_pathways_up_per_module.csv
        │   │   ├── Enrichr_top_10_enriched_pathways_down_per_module.csv
        │   │   ├── module_1/
        │   │   │   ├── GO_Biological_Process_2025.png
        │   │   │   ├── ...
        │   │   └── ...
        │   │
        │   ├── TFs2targets/
        │   │   ├── TP53/
        │   │   │   └── GO_Biological_Process_2025.png
        │   │   └── ...
        │   ├── Metabolites2targets/
        │   │   ├── L_glutamine/
        │   │   │   └── GO_Biological_Process_2025.png
        │   │   └── ...
        │   │
        │   └── Other_omics2targets/
        │       ├── feature/
        │       │   └── GO_Biological_Process_2025.png
        │       └── ...
        │
        └── Module_Overview/
            ├── Module_Overview.csv (main result summary file)
            ├── interactive_module_network.html
            ├── Module_Expression_Heatmap.png
            ├── module_expression_analysis.csv
            ├── megago_similarity_matrix.csv
            └── megaGO_files/ (used for megago clustering)
```

### Key Output File Descriptions

#### Preprocessing Files

**`LemonPreprocessed_expression_TFA.txt`**
- Normalized gene expression + TFA scores
- Format: Genes/TFs (rows) × Samples (columns)
- Used for LemonTree clustering

**`LemonPreprocessed_metabolomics.txt`**
- Normalized and scaled metabolite abundances
- Pareto scaling applied (if `use_omics_specific_scaling=true`)

**`LemonPreprocessed_complete.txt`**
- Combined matrix: genes + TFA + metabolites + other omics
- Input to LemonTree clustering

**`TFA_consensus.txt`**
- Transcription factor activity scores
- Generated by decoupleR
- One row per TF, one column per sample

#### Network Files

**`LemonNetwork_{method}_{N}modules.txt`**
- Complete regulatory network
- Columns:
  - `Regulator`: Regulator name
  - `Target`: Gene name
  - `Score`: LemonTree assignment score
  - `Lemon_module`: Module ID
  - `Type`: Regulator type (TF, Metabolite, etc.)

**`{Prefix}2targets_{method}_{N}_modules.txt`**
- Per-regulator-type edge lists
- Columns:
  - `Regulator`: Name
  - `Targets`: Pipe-separated list of target genes
- Example filename: `Lovering2targets_fold2.5x_permodule_15_modules.txt`

**Cytoscape Files:**
- `Cytoscape_edges.txt`: All edges for network import into cytoscape
- `Cytoscape_nodes.txt`: Node attributes (type, module, etc.) to create nice visualizations
- Per-module files in `subnetworks/cytoscape/` to visualize individual modules in the Lemonite knowledge graph

#### Subnetwork Graphs

**`graph_{module}_categorized_clean.png`**
- Per-module network visualization
- Node colors:
  - Red: Metabolites/Lipids
  - Green: Transcription factors
  - Orange: Target genes
- Edge types:
  - Arrows: Causal (LINCS, chEMBL)
  - Circles: Metabolic pathway (GEM)
  - Lines: PPI (STITCH, BioGRID, IntAct, STRING)
  - Dashed: Ambiguous

#### Module Overview

**`Module_Overview.csv`**
Master table with columns:
- `Module`: Module ID
- `Module_genes`: Module genes separeted by |
- `Coherence`: Module coherence score
- `TF_regulators`: TF regulators for module
- `Metabolite_regulators`: Metabolite regulators for module
- `Other_omics_regulators`: Other omics regulators, 1 column per other omics provided
- `Top_3_pathways_molecular_function`: Top 3 enriched pathways in GO Molecular Function 2025
- `Top_3_pathways_cellular_component`: Top 3 enriched pathways in GO Cellular Component 2025
- `Top_3_pathways_biological_process`: Top 3 enriched pathways in GO Biological Process 2025
- `Top_3_pathways_Reactome`: Top 3 enriched pathways in Reactome 2025
- `Top_3_pathways_KEGG`: Top 3 enriched pathways in KEGG 2021
- `Expression_adjustedpval`: adjusted p-value for module differential expression

**`interactive_module_network.html`**
- Interactive Plotly-based visualization of the integrated regulator-module network
- Additional information displayed when hovering over a node

---

## Advanced Usage

### Custom Regulator Types

#### Adding New Regulator Types

The pipeline supports any number and type of regulators through the `--regulator_types` parameter:

```bash
# Example: Add proteins as regulators
--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Proteins:Proteomics.txt"

# Example: Add kinases and microRNAs + use a custom TF list
--regulator_types "TFs:MyTFs.txt,Kinases:Kinase_abundance.txt,miRNAs:miRNA_expression.txt"
```

**Requirements:**
1. Create your data file in the `data/` directory with the exact name specified
2. For abundance data: Format as tab-separated with features as rows, samples as columns
3. For TF lists: Format as plain text with one gene symbol per line

**TF List File Example (`data/MyTFs.txt`):**
```txt
TP53
MYC
STAT3
NFKB1
```

**Abundance File Example (`data/Proteomics.txt`):**
```tsv
Name	Sample1	Sample2	Sample3
AKT1	234.5	345.6	456.7
MAPK1	567.8	678.9	789.0
CDK1	123.4	234.5	345.6
```

**How it works:**
1. **Preprocessing stage**: 
   - For TF lists: Copied to `Preprocessing/` directory with lowercase prefix (e.g., `tfs.txt`)
   - For abundance data: Normalized, scaled, and row names extracted to create regulator list (e.g., `proteins.txt`)
   - All abundance data merged into `LemonPreprocessed_complete.txt`

2. **LemonTree clustering**: 
   - Runs regulator assignment for each regulator type
   - Creates `{Prefix}.allreg.txt` and `{Prefix}.randomreg.txt` files

3. **Network generation**: 
   - Selects top regulators per module
   - Creates `{Prefix}2targets_*.txt` edge lists
   - Combines all into unified network with type annotations

### Custom TFA Network

By default we use the curated [CollecTri network](https://academic.oup.com/nar/article/51/20/10934/7318114) of TF-target gene interactions for TFA inference. You can also provide your own network:

```bash
--prior_network /path/to/custom_network.tsv
```

**Format:** Tab-separated
- Column 1: `source` (TF gene symbol)
- Column 2: `target` (target gene symbol)
- Optional Column 3: `weight` (interaction confidence)

**Example:**
```tsv
source	target
TP53	CDKN1A
TP53	MDM2
MYC	CCND1	
```


### Fine-Tuning LemonTree Parameters

```bash
# High-quality analysis (slower)
--n_clusters 200 \
--coherence_threshold 0.7 \
--regulator_fold_cutoff 3.0

# Fast exploratory analysis
--n_clusters 20 \ # not recommened! 
--coherence_threshold 0.4 \
--regulator_selection_method percentage \
--top_n_percent_regulators 5.0

# Strict filtering
--coherence_threshold 0.8 \
--regulator_fold_cutoff 4.0 \
--top_n_genes 1000
```

### Working with Pre-normalized Data

If your omics data is already normalized and log-transformed, the pipeline should auto detect this (median < 100) and skip log-transformation.

---

### Module Heatmap Annotations

The pipeline can display multiple metadata annotations as colored bars on module heatmaps, allowing you to visualize sample groupings by different clinical or experimental variables.

#### Overview

**What metadata is available?**
- All columns saved during preprocessing via `--metadata_columns` parameter are stored in:
  - `LemonTree/Preprocessing/DESeq_groups.txt` (used by DESeq2)
  - `ModuleViewer_files/sample_mapping.mvf` (used by module_viewer.py)

**How to control what's displayed?**
- By default, only the primary contrast variable (`diagnosis`) is shown as an annotation bar
- Use the module_viewer.py `--annotation_types` parameter to select which annotations to display
- This is **independent** of `--metadata_columns` and `--design_formula` - you can save many columns but display only a subset

#### Configuration

**Step 1: Save Metadata During Preprocessing**

Specify which metadata columns to include when running the pipeline:

```bash
nextflow run main.nf \
  --input_dir /path/to/data \
  --metadata_columns "diagnosis,sex,age,batch,tissue" \
  -profile docker
```

This saves all specified columns in `DESeq_groups.txt` and generates `sample_mapping.mvf` with color mappings for each metadata type.

**Step 2: Select Annotations for Heatmaps**

The module viewer automatically displays annotation bars based on available metadata. To customize which annotations are shown, you can modify the viewer_heatmaps module parameters (see Advanced Customization below).

By default, the pipeline will:
- Display the `diagnosis` annotation bar (primary contrast variable)
- Show all other available metadata as additional annotation bars

#### Annotation Bar Display

**Example with Multiple Annotations:**

```
┌─────────────────────────────────────────────┐
│   TFs Heatmap                               │  ← Regulator heatmap
├─────────────────────────────────────────────┤
│   Metabolites Heatmap                       │  ← Regulator heatmap
├─────────────────────────────────────────────┤
│   Expression Data Heatmap                   │  ← Gene expression heatmap
├─────────────────────────────────────────────┤
│ [RED][RED][BLUE][BLUE][RED]...             │  ← Diagnosis (e.g., case/control)
├─────────────────────────────────────────────┤
│ [BLUE][PINK][BLUE][PINK][BLUE]...          │  ← Sex (e.g., M/F)
├─────────────────────────────────────────────┤
│ [GREEN][GREEN][ORANGE][ORANGE][GREEN]...   │  ← Batch (e.g., 1/2/3)
├─────────────────────────────────────────────┤
│ Legend:                                     │  ← Combined legend
│  Diagnosis: RED=case  BLUE=control          │
│  Sex: BLUE=Male  PINK=Female                │
│  Batch: GREEN=1  ORANGE=2  PURPLE=3         │
└─────────────────────────────────────────────┘
```

#### Advanced Customization

**Custom Annotation Selection (Direct module_viewer.py call):**

If running module_viewer.py separately (not through Nextflow), you can specify exactly which annotations to display:

```bash
python scripts/module_viewer.py \
  --input_dir results/my_run/LemonTree \
  --output_dir results/my_run/custom_heatmaps \
  --regulator_files "TFs:TFs.selected_regs_list.txt,Metabolites:Metabolites.selected_regs_list.txt" \
  --annotation_types "diagnosis,sex,batch" \
  --annotation_labels "Disease Status,Biological Sex,Experimental Batch"
```

**Parameters:**
- `--annotation_types`: Comma-separated list of metadata column names to display (must be in sample_mapping.mvf)
- `--annotation_labels`: Optional custom display labels (must match number of annotation types)

**Examples:**

```bash
# Display only diagnosis (default behavior)
--annotation_types "diagnosis"

# Display diagnosis and sex
--annotation_types "diagnosis,sex"

# Display multiple annotations with custom labels
--annotation_types "diagnosis,sex,batch,tissue" \
--annotation_labels "Clinical Status,Biological Sex,Batch ID,Tissue Type"

# Display only batch (useful for batch effect visualization)
--annotation_types "batch"
```

#### Color Schemes

The pipeline automatically assigns colors to metadata values:

**Categorical Variables:**
- General categories: RED, BLUE, GREEN, ORANGE, PURPLE, YELLOW, PINK, CYAN, etc.
- Cycles through palette if more than 16 categories

**Sex/Gender (Special Handling):**
- Male/M: BLUE
- Female/F: PINK
- Unknown: GREY

**Automatic Detection:**
- Columns with >20 unique values are skipped (assumed continuous)
- Columns with all NA values are skipped
- Numeric columns are treated as categorical (create bins or skip if too many values)

#### File Formats

**sample_mapping.mvf Format (Multi-Metadata):**

The MVF file contains multiple sections separated by `---`, one for each metadata type:

```
::TYPE=diagnosis
::GLOBAL
::VALUES=color
::OBJECT=CONDITIONS
::LEGEND=case:RED|control:BLUE	Clinical Status
|Sample1:RED|Sample2:BLUE|Sample3:RED|Sample4:BLUE
---
::TYPE=sex
::GLOBAL
::VALUES=color
::OBJECT=CONDITIONS
::LEGEND=M:BLUE|F:PINK	Biological Sex
|Sample1:BLUE|Sample2:PINK|Sample3:BLUE|Sample4:PINK
---
::TYPE=batch
::GLOBAL
::VALUES=color
::OBJECT=CONDITIONS
::LEGEND=1:GREEN|2:ORANGE|3:PURPLE	Batch
|Sample1:GREEN|Sample2:GREEN|Sample3:ORANGE|Sample4:PURPLE
```

**Backward Compatibility:**
- Old single-metadata MVF files (without `---` separator) are still supported
- Automatically detected and parsed correctly

#### Troubleshooting

**"Requested annotation type 'X' not found":**
- Check that the column was included in `--metadata_columns` during preprocessing
- Verify the column name matches exactly (case-sensitive)
- Check `LemonTree/Preprocessing/DESeq_groups.txt` to see available columns

**Too many unique values warning:**
- Columns with >20 unique values are automatically skipped
- For continuous variables (age, score), consider binning them in your metadata file
- Example: Replace `age` column (25, 30, 35...) with `age_group` (<30, 30-50, >50)

**Colors not displaying correctly:**
- Check `ModuleViewer_files/sample_mapping.mvf` to verify color assignments
- Colors are normalized to lowercase (RED → red, BLUE → blue)
- Unknown/NA values default to grey

**No annotation bars displayed:**
- Ensure `sample_mapping.mvf` exists in `ModuleViewer_files/`
- Check pipeline logs for preprocessing warnings
- Verify at least one metadata column has valid (non-NA) values

#### Best Practices

1. **Include relevant metadata early**: Specify `--metadata_columns` with all potentially useful columns during pipeline run
2. **Keep categories manageable**: Limit categorical variables to <20 unique values
3. **Use meaningful names**: Column names like `tissue` or `biopsy_location` are more interpretable than `var1` or `group_x`
4. **Bin continuous variables**: Convert age, BMI, or scores into categorical groups before analysis
5. **Custom labels for clarity**: Use `--annotation_labels` to provide publication-ready annotation names

#### Example Workflow

```bash
# 1. Run pipeline with multiple metadata columns
nextflow run main.nf \
  --input_dir /data/sepsis_study \
  --metadata_columns "diagnosis,sex,age_group,batch,tissue" \
  --design_formula "~ batch + diagnosis" \
  -profile docker

# Pipeline automatically generates module heatmaps with all available annotations

# 2. (Optional) Generate custom heatmaps with selected annotations
python scripts/module_viewer.py \
  --input_dir results/sepsis_study/LemonTree \
  --output_dir results/sepsis_study/custom_viz \
  --regulator_files "TFs:TFs.selected_regs_list.txt,Metabolites:Metabolites.selected_regs_list.txt" \
  --annotation_types "diagnosis,tissue" \
  --annotation_labels "Sepsis Status,Tissue Origin"
```

---

## Performance Optimization

### Parallelization Strategies

#### Local/Docker
```bash
# Optimize for system with 16 cores, 64GB RAM
--max_cpus 16 \
--max_memory 64.GB \
-profile docker

# In nextflow.config or via CLI
process.executor = 'local'
executor.queueSize = 50  # Allow up to 50 parallel tasks
```

#### HPC/Singularity
```bash
# Standard HPC with Singularity
--max_cpus 20 \
--max_memory 128.GB \
-profile singularity

# Optimized HPC profile (recommended for large datasets)
-profile hpc,singularity \
--max_cpus 32 \
--max_memory 128.GB

# In nextflow.config (for SLURM):
process.executor = 'slurm'
process.queue = 'normal'
process.clusterOptions = '--account=myproject'
```

**HPC Profile Features:**
- Extended time limits (up to 72h for intensive processes)
- Increased memory and CPU allocations
- Aggressive retry strategy for common HPC exit codes (143, 137, etc.)
- Optimized for large-scale analyses with megago clustering

**Using MegaGO on HPC:**
```bash
# Enable megago clustering (computationally intensive!)
nextflow run main.nf \
  -profile hpc,singularity \
  --input_dir ./your_data \
  --use_megago true \
  --overview_n_clusters 5 \
  [other parameters]
```

### Resource Allocation Guidelines

| Dataset Size | Recommended Resources |
|--------------|----------------------|
| Small (<50 samples, <10K genes) | 4 cores, 16GB RAM |
| Medium (50-200 samples, 10-20K genes) | 8 cores, 32GB RAM |
| Large (>200 samples, >20K genes) | 16+ cores, 64GB+ RAM |

---

## Container Management

### Docker Image Management

#### Building
```bash
# Standard build
docker build -t lemontree-pipeline:latest .

# Build without cache
docker build --no-cache -t lemontree-pipeline:latest .

# Build with specific tag
docker build -t lemontree-pipeline:v1.0.0 .
```

#### Listing
```bash
# List all images
docker images

# List lemontree images
docker images | grep lemontree-pipeline
```

#### Removing
```bash
# Remove specific image
docker rmi lemontree-pipeline:latest

# Force remove (if containers exist)
docker rmi -f lemontree-pipeline:latest

# Remove all lemontree images
docker rmi $(docker images | grep lemontree-pipeline | awk '{print $3}')
```

#### Cleaning Up
```bash
# Stop all containers
docker stop $(docker ps -a -q)

# Remove all stopped containers
docker container prune

# Remove unused images
docker image prune -a

# Complete cleanup (WARNING: removes all unused Docker data)
docker system prune -a
```

### Singularity Image Management

#### Building
```bash
# Standard build (requires sudo)
sudo singularity build lemontree-pipeline.sif Singularity.def

# Build in sandbox mode (for testing)
sudo singularity build --sandbox lemontree-sandbox/ Singularity.def

# Build from Docker image
singularity build lemontree-pipeline.sif docker://lemontree-pipeline:latest
```

#### Testing

Run the script ./test_singularity.sh to build the singularity image and do a Lemonite run on the test dataset provided in this repository.


---

## FAQ

### General Questions

**Q: How long does the pipeline take?** 
A: This mostly depends on the number of genes included (HVG selection threshold) and the number of LemonTree clustering runs.

**Q: Can I run without metabolomics data?**
A: Yes! The pipeline is flexible. You can use any combination of regulator types. For TFs only: `--regulator_types "TFs:Lovering_TF_list.txt"`. Just know that the evaluation process against the Lemonite knowledge graph was developed for metabolites.

**Q: What if my metabolomics data uses different IDs (KEGG, ChEBI)?**
A: Convert to HMDB IDs using MetaboAnalyst or create custom mapping CSV. Note: This is only needed if you want to evaluate against the Lemonite knowledge graph.

### Technical Questions

**Q: Why use Pareto scaling for metabolomics?**
A: Pareto scaling preserves variance structure better than z-score, maintains relative magnitude differences, and is less sensitive to outliers—all important for metabolomics data.

**Q: Can I use proteomics data?**
A: Yes! Just provide the file and specify in regulator_types: `--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Proteins:Proteomics.txt"`

**Q: How do I choose coherence_threshold?**
A: Start with 0.6 (default). Check `Module_coherence_scores.txt` to see distribution. Lower if too few modules, raise if too many low-quality modules.

**Q: I ran Lemonite twice and the results are different?**
A: The LemonTree algorithm does have some stochasticity, so that differences between individual LemonTree runs can occur.

**Q: Pipeline fails with "IndexOutOfBoundsException" during POST_CLUSTERING for mouse data?**
A: This error occurs when the pipeline uses human TF symbols with mouse expression data, resulting in zero TF matches. **Solution:**
1. Ensure you specify `--organism mouse` when running the pipeline
2. The pipeline will automatically use `lovering_TF_list_mouse.txt` (1,141 mouse TF symbols)
3. If you get this error, the mouse TF list should already be in `PKN/lovering_TF_list_mouse.txt` and will be automatically copied to your data directory
4. Rerun with `-resume` to continue from where it failed: `nextflow run main.nf --input_dir /path/to/data --organism mouse -profile singularity -resume`

**Q: How do I know if I should use --organism mouse or --organism human?**
A: Check your expression data file (`Counts.tsv`). If gene symbols look like `0610009L18Rik`, `Slc25a5`, `Gm42418` (mixed case, often with numbers), you have mouse data. If they look like `TP53`, `EGFR`, `MYC` (all uppercase), you have human data.

---

## Citation & References

### Lemonite Citation
```
Vandemoortele et al. (2025). Lemonite: Interpretable multi-omics integration for
regulatory metabolite discovery. [Journal], [Volume], [Pages].
```

### LemonTree Algorithm
```
Bonnet, E., Calzone, L., & Michoel, T. (2015). Integrative Multi-omics Module Network
Inference with Lemon-Tree. PLOS Computational Biology, 11(2), e1003983.
https://doi.org/10.1371/journal.pcbi.1003983
```

---

## Contact & Support

- **Issues**: [GitHub Issues](https://github.com/CBIGR/Lemonite/issues)
- **Email**: boris.vandemoortele@ugent.be
- **Lab**: CBIGR lab @ Ghent University

---

**Last Updated:** November 26, 2025
