#!/bin/bash
#
# Lemonite Pipeline - Sepsis Dataset Analysis
# Dataset: /home/borisvdm/Documents/PhD/25_sepsis_Libert
# 
# Configuration:
# - 10 clusters
# - 5 CPUs
# - 3000 highly variable genes
# - Regulators: TFs, metabolites, proteins
# - Metadata annotations: diagnosis, survival, sex, smoker, alcoholic

set -e  # Exit on error
set -u  # Exit on undefined variable

# ============================================================================
# CONFIGURATION
# ============================================================================

INPUT_DIR="/home/borisvdm/Documents/PhD/25_sepsis_Libert"
PIPELINE_DIR="/home/borisvdm/repo/LemonIte/nextflow"

# Pipeline parameters
N_CLUSTERS=10
MAX_CPUS=10
TOP_N_GENES=3000

# Regulator types with their data files
# Format: "Prefix:DataFile,Prefix:DataFile"
# Note: File names are case-sensitive and must match exactly
REGULATOR_TYPES="TFs:lovering_TF_list.txt,Metabolites:Metabolomics.txt,Proteins:proteomics.txt"

# Metadata columns for heatmap annotations
HEATMAP_METADATA="diagnosis,survival,sex,smoker,alcoholic"

# DESeq2 and differential expression settings
DESEQ_CONTRAST="diagnosis"
DESIGN_FORMULA="~ diagnosis"
METADATA_COLS="diagnosis,survival,sex,smoker,alcoholic"

# Optional: Specify run identifier (leave empty for auto-generation)
RUN_ID=""  # Will be auto-generated based on parameters

# ============================================================================
# PRE-FLIGHT CHECKS
# ============================================================================

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                                                                ║"
echo "║     Lemonite Pipeline Setup - Sepsis Dataset                   ║"
echo "║                                                                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "❌ ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Check if data subdirectory exists
if [ ! -d "$INPUT_DIR/data" ]; then
    echo "❌ ERROR: Data subdirectory not found: $INPUT_DIR/data"
    echo "   Please ensure your input data is organized in $INPUT_DIR/data/"
    exit 1
fi

# Check if pipeline directory exists
if [ ! -d "$PIPELINE_DIR" ]; then
    echo "❌ ERROR: Pipeline directory not found: $PIPELINE_DIR"
    exit 1
fi

# Check if main.nf exists
if [ ! -f "$PIPELINE_DIR/main.nf" ]; then
    echo "❌ ERROR: main.nf not found in $PIPELINE_DIR"
    exit 1
fi

echo "✅ All pre-flight checks passed!"
echo ""

# ============================================================================
# DISPLAY CONFIGURATION
# ============================================================================

echo "📋 Pipeline Configuration:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Input Directory:        $INPUT_DIR"
echo "  Data Directory:         $INPUT_DIR/data"
echo "  Output Directory:       $INPUT_DIR/results"
echo "  Work Directory:         $INPUT_DIR/work"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Number of Clusters:     $N_CLUSTERS"
echo "  Maximum CPUs:           $MAX_CPUS"
echo "  Top Variable Genes:     $TOP_N_GENES"
echo "  Regulator Types:        $REGULATOR_TYPES"
echo "  Heatmap Metadata:       $HEATMAP_METADATA"
echo "  DESeq Contrast:         $DESEQ_CONTRAST"
echo "  Design Formula:         $DESIGN_FORMULA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# ============================================================================
# REQUIRED INPUT FILES CHECK
# ============================================================================

echo "🔍 Checking required input files..."

REQUIRED_FILES=(
    "data/Counts.tsv"
    "data/Metadata.txt"
)

OPTIONAL_FILES=(
    "data/lovering_TF_list.txt"
    "data/Metabolomics.txt"
    "data/proteomics.txt"
    "data/Lipidomics.txt"
    "data/name_map.csv"
)

ALL_OK=true

for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$INPUT_DIR/$file" ]; then
        echo "  ✅ Found: $file"
    else
        echo "  ❌ Missing: $file"
        ALL_OK=false
    fi
done

echo ""
echo "Optional files (specified as regulators):"
for file in "${OPTIONAL_FILES[@]}"; do
    if [ -f "$INPUT_DIR/$file" ]; then
        echo "  ✅ Found: $file"
    else
        echo "  ⚠️  Not found: $file (will be skipped in analysis)"
    fi
done
echo ""

if [ "$ALL_OK" = false ]; then
    echo "❌ ERROR: Required input files are missing!"
    echo ""
    echo "Required file structure:"
    echo "$INPUT_DIR/"
    echo "├── data/"
    echo "│   ├── Counts.tsv               # REQUIRED: Gene expression matrix (raw counts)"
    echo "│   ├── Metadata.txt             # REQUIRED: Sample metadata (must include: $METADATA_COLS)"
    echo "│   ├── Lovering_TF_list.txt     # OPTIONAL: List of transcription factors"
    echo "│   ├── Metabolomics.txt         # OPTIONAL: Metabolomics data"
    echo "│   ├── Proteomics.txt           # OPTIONAL: Proteomics data"
    echo "│   ├── Lipidomics.txt           # OPTIONAL: Lipidomics data"
    echo "│   └── name_map.csv             # OPTIONAL: Metabolite ID mapping (HMDB)"
    echo ""
    echo "Note: File names are case-sensitive!"
    echo ""
    exit 1
fi

# ============================================================================
# RUN PIPELINE
# ============================================================================

echo "🚀 Starting Lemonite pipeline..."
echo ""

cd "$PIPELINE_DIR"

# Build the nextflow command
NEXTFLOW_CMD="nextflow run main.nf \
  --input_dir $INPUT_DIR \
  --n_clusters $N_CLUSTERS \
  --max_cpus $MAX_CPUS \
  --top_n_genes $TOP_N_GENES \
  --regulator_types \"$REGULATOR_TYPES\" \
  --heatmap_metadata_cols \"$HEATMAP_METADATA\" \
  --deseq_contrast1 \"$DESEQ_CONTRAST\" \
  --design_formula \"$DESIGN_FORMULA\" \
  --metadata_columns \"$METADATA_COLS\" \
  --perform_tfa true \
  --interactive_overview true \
  --use_megago false \
  -profile singularity \
  -resume"

# Add run_id if specified
if [ -n "$RUN_ID" ]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD --run_id $RUN_ID"
fi

echo "Executing command:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "$NEXTFLOW_CMD"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Execute the pipeline
eval $NEXTFLOW_CMD

# ============================================================================
# POST-RUN INFORMATION
# ============================================================================

if [ $? -eq 0 ]; then
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                                                                ║"
    echo "║          ✅ Pipeline completed successfully! ✅                ║"
    echo "║                                                                ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "📁 Results are available in: $INPUT_DIR/results"
    echo ""
    echo "📊 Key output files:"
    echo "  • Module overview:     results/Module_overview/module_overview.html"
    echo "  • Network files:       results/Networks/"
    echo "  • Enrichment results:  results/Enrichment/"
    echo "  • Module heatmaps:     results/Module_viewer/"
    echo "  • PKN evaluation:      results/PKN_evaluation/"
    echo ""
    echo "🔬 Analysis features enabled:"
    echo "  ✓ Transcription Factor Activity (TFA) calculation"
    echo "  ✓ Multi-regulator integration (TFs, metabolites, proteins)"
    echo "  ✓ Interactive module overview with functional clustering"
    echo "  ✓ Differential expression analysis ($DESEQ_CONTRAST)"
    echo ""
    echo "ℹ️  Note: MegaGO clustering is disabled by default (use --use_megago true to enable)"
    echo ""
else
    echo ""
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                                                                ║"
    echo "║          ❌ Pipeline failed! ❌                                ║"
    echo "║                                                                ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "🔍 Check the logs:"
    echo "  • Nextflow log: .nextflow.log"
    echo "  • Work directory: $INPUT_DIR/work"
    echo ""
fi
