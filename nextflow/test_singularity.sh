#!/bin/bash

#############################################################################
# Lemonite Pipeline Test Script (Singularity)
#############################################################################
# This script builds the Singularity image and runs the Lemonite pipeline
# with test data for validation purposes.
#
# Author: Vandemoortele Boris, CBIGR lab @ Ghent University
# Date: November 18, 2025
#############################################################################

set -e  # Exit on any error

echo "╔════════════════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                                ║"
echo "║                L E M O N I T E   P I P E L I N E   T E S T   ( S I N G U L A R I T Y )         ║"
echo "║                                                                                                ║"
echo "╚════════════════════════════════════════════════════════════════════════════════════════════════╝"
echo ""

#############################################################################
# Configuration
#############################################################################
PIPELINE_DIR=${PWD}
INPUT_DIR=${PWD}"/test_dataset"
SINGULARITY_IMAGE="lemontree-pipeline.sif"
N_CLUSTERS=8
TOP_N_GENES=1500
MAX_CPUS=8

echo "Configuration:"
echo "   Pipeline directory: $PIPELINE_DIR"
echo "   Input directory: $INPUT_DIR"
echo "   Singularity image: $SINGULARITY_IMAGE"
echo "   Number of clusters: $N_CLUSTERS"
echo "   Top N genes: $TOP_N_GENES"
echo "   Max CPUs: $MAX_CPUS"
echo ""

#############################################################################
# Step 1: Check prerequisites
#############################################################################
echo "Step 1: Checking prerequisites..."

# Check if Singularity is installed
if ! command -v singularity &> /dev/null; then
    echo "ERROR: Singularity is not installed or not in PATH"
    echo "   Please install Singularity from: https://sylabs.io/guides/latest/user-guide/"
    exit 1
fi
echo "   Singularity found: $(singularity --version)"

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: Nextflow is not installed or not in PATH"
    echo "   Please install Nextflow from: https://www.nextflow.io/"
    exit 1
fi
echo "   Nextflow found: $(nextflow -version | head -1)"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi
echo "   Input directory found"

# Check if pipeline directory exists
if [ ! -d "$PIPELINE_DIR" ]; then
    echo "ERROR: Pipeline directory not found: $PIPELINE_DIR"
    exit 1
fi
echo "   Pipeline directory found"

echo ""

#############################################################################
# Step 2: Build Singularity image (if not exists)
#############################################################################
echo "Step 2: Building Singularity image..."

cd "$PIPELINE_DIR"

if [ -f "$SINGULARITY_IMAGE" ]; then
    echo "   INFO: Singularity image already exists: $SINGULARITY_IMAGE"
    echo "   Using existing image"
else
    echo "   Building Singularity image for the first time..."
    echo "   WARNING: This may take about 2 hours..."
    ./build-singularity.sh
fi

if [ ! -f "$SINGULARITY_IMAGE" ]; then
    echo "Error: Failed to build Singularity image"
    exit 1
fi
echo "   Singularity image ready: $SINGULARITY_IMAGE"

echo ""

#############################################################################
# Step 3: Run the pipeline
#############################################################################
echo "Step 3: Running Lemonite pipeline..."
echo ""
echo "   This will run with the following parameters:"
echo "   - Input: $INPUT_DIR"
echo "   - Regulator types: TFs, Metabolites, Lipids, Proteins"
echo "   - Number of clusters: $N_CLUSTERS"
echo "   - Top N genes: $TOP_N_GENES"
echo "   - Enrichment method: enrichr"
echo "   - Max CPUs: $MAX_CPUS"
echo ""
echo "   Estimated time: about 2 hours for building singularity image, 30-60 minutes for running pipeline (depending on system)"
echo ""

# Run the pipeline
nextflow run "$PIPELINE_DIR/main.nf" \
  --input_dir "$INPUT_DIR" \
  --regulator_types "TFs:lovering_TF_list.txt,Metabolites:Metabolomics.txt,Lipids:lipidomics.tsv,Proteins:Proteomics_TFs.txt" \
  --enrichment_method "enrichr" \
  --n_clusters "$N_CLUSTERS" \
  --top_n_genes "$TOP_N_GENES" \
  -profile singularity \
  --max_cpus "$MAX_CPUS" \
  -resume

#############################################################################
# Step 4: Check results
#############################################################################
echo ""
echo "╔════════════════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                                ║"
echo "║                             P I P E L I N E   C O M P L E T E D                                ║"
echo "║                                                                                                ║"
echo "╚════════════════════════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Results summary:"
echo ""

# Check if results directory exists
RESULTS_DIR="${INPUT_DIR}/results"
if [ -d "$RESULTS_DIR" ]; then
    echo "   Results directory: $(realpath "$RESULTS_DIR")"
    echo ""
    echo "   Output files:"
    
    # Find the run directory (should match the pattern *HVG_coherence*)
    RUN_DIR=$(find "$RESULTS_DIR" -maxdepth 1 -type d -name "*HVG_coherence*" | head -1)
    
    if [ -n "$RUN_DIR" ]; then
        echo "   Run directory: $(basename "$RUN_DIR")"
        echo ""
        
        # List key output files
        if [ -d "$RUN_DIR/LemonTree" ]; then
            echo "      - LemonTree files: $(find "$RUN_DIR/LemonTree" -type f 2>/dev/null | wc -l) files"
        fi
        
        if [ -d "$RUN_DIR/Networks" ]; then
            echo "      - Network files: $(find "$RUN_DIR/Networks" -type f -name "*.tsv" 2>/dev/null | wc -l) files"
        fi
        
        if [ -d "$RUN_DIR/module_heatmaps" ]; then
            echo "      - Module heatmaps: $(find "$RUN_DIR/module_heatmaps" -type f -name "*.pdf" 2>/dev/null | wc -l) files"
        fi
        
        if [ -d "$RUN_DIR/enrichment" ]; then
            echo "      - Enrichment results: $(find "$RUN_DIR/enrichment" -type f -name "*.csv" 2>/dev/null | wc -l) files"
        fi
        
        if [ -f "$RUN_DIR/module_overview/interactive_module_overview.html" ]; then
            echo "      - Interactive overview: $(basename "$RUN_DIR")/module_overview/interactive_module_overview.html"
            echo ""
            echo "   Open the interactive overview in your browser:"
            echo "      file://$(realpath "$RUN_DIR/module_overview/interactive_module_overview.html")"
        fi
    fi
else
    echo "   WARNING: Results directory not found at $RESULTS_DIR"
fi

echo ""
echo "Test completed successfully!"
echo ""
