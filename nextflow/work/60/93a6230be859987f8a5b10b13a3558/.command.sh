#!/bin/bash -euo pipefail
# Create working directory
workdir=$(pwd)

# Find required files
cluster_file=$(find . -name "clusters_list.txt" | head -1)
coherence_file="Module_coherence_scores.txt"

# Use the provided expression and DESeq files
expression_file="LemonPreprocessed_expression.txt"
deseq_groups_file="DESeq_groups.txt"

# Also look for the MVF file that should be created by PKN evaluation
mvf_file=$(find . -name "metabolite_LemoniteKG_interactions.mvf" | head -1)

echo "=== Module Overview Generation ==="
echo "Found files:"
echo "  Clusters: $cluster_file"
echo "  Coherence scores: $coherence_file"
echo "Expression: $expression_file (provided as input)"
echo "  DESeq groups: $deseq_groups_file (provided as input)"
echo "  MVF file: $mvf_file"
echo "Regulator types: TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"
echo "Parameters:"
echo "  Enrichment method: EnrichR"
echo "  N clusters: 5"
echo "  Organism: human"
echo "  Coherence threshold: 0.5"
echo "  Clustering workflow: canonical MegaGO top_30 + rrvgo"
echo "  PKN file: Lemonite_PKN.tsv"
echo "  Metabolite mapping: name_map.csv"
# Debug: list enrichment dir if present
if [ -d "Enrichment" ]; then
    echo "=== Enrichment directory contents ==="
    ls -la Enrichment || true
    echo "====================================="
fi

# Build regulator_files argument dynamically based on regulator_types
REGULATOR_FILES=""
REGULATOR_SCORE_FILES=""
REGULATOR_TYPES="TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"

# Parse regulator types and find corresponding files
IFS=',' read -ra TYPE_ARRAY <<< "$REGULATOR_TYPES"
for type_config in "${TYPE_ARRAY[@]}"; do
    # Extract type name (before colon) - this is the Prefix like "TFs", "Metabolites"
    TYPE_NAME=$(echo "$type_config" | cut -d':' -f1)
    # Extract data filename (after colon) - not used for file lookup
    FILE_PREFIX=$(echo "$type_config" | cut -d':' -f2)

    # Find the corresponding regulator file using TYPE_NAME (Prefix), NOT data filename
    # Files are named like: TFs.selected_regs_list.txt, Metabolites.selected_regs_list.txt
    REG_FILE="${TYPE_NAME}.selected_regs_list.txt"
    if [ -f "$REG_FILE" ]; then
        if [ -n "$REGULATOR_FILES" ]; then
            REGULATOR_FILES="$REGULATOR_FILES,${TYPE_NAME}:${REG_FILE}"
        else
            REGULATOR_FILES="${TYPE_NAME}:${REG_FILE}"
        fi
        echo "Found ${TYPE_NAME} regulator file: $REG_FILE"
    else
        echo "Warning: Regulator file $REG_FILE not found"
    fi

    # Also look for regulator score files (selected_regulators_scores.txt)
    SCORE_FILE="${TYPE_NAME}.selected_regulators_scores.txt"
    if [ -f "$SCORE_FILE" ]; then
        if [ -n "$REGULATOR_SCORE_FILES" ]; then
            REGULATOR_SCORE_FILES="$REGULATOR_SCORE_FILES,${TYPE_NAME}:${SCORE_FILE}"
        else
            REGULATOR_SCORE_FILES="${TYPE_NAME}:${SCORE_FILE}"
        fi
        echo "Found ${TYPE_NAME} regulator score file: $SCORE_FILE"
    else
        echo "Warning: Regulator score file $SCORE_FILE not found"
    fi
done

# Construct command arguments
args="--input_dir . --output_dir ."
args="$args --enrichment_method EnrichR"
args="$args --coherence_threshold 0.5"
args="$args --group_column diagnosis"
if [ -n "$REGULATOR_FILES" ]; then
    args="$args --regulator_files $REGULATOR_FILES"
fi
if [ -n "$REGULATOR_SCORE_FILES" ]; then
    args="$args --regulator_score_files $REGULATOR_SCORE_FILES"
fi

# Always pass n_clusters for canonical MegaGO top_30 clustering
args="$args --n_clusters 5"
args="$args --organism human"
echo "Functional clustering: canonical MegaGO top_30, n_clusters=5"

# Add PKN file for edge categorization if available
if [ -f "Lemonite_PKN.tsv" ]; then
    args="$args --pkn_file Lemonite_PKN.tsv"
    echo "PKN file provided for edge categorization"
fi

# Add metabolite mapping if available
if [ -f "name_map.csv" ]; then
    args="$args --metabolite_mapping name_map.csv"
    echo "Metabolite name mapping provided"
fi

# Add name mapping for restoring original feature names
if [ -f "name_mapping.tsv" ]; then
    args="$args --name_mapping name_mapping.tsv"
    echo "Name mapping provided for original name restoration"
fi

# Expression prioritization is enabled by default but can be disabled
if [ "true" = "true" ]; then
    args="$args --prioritize_by_expression"
    args="$args --expression_file $expression_file"
    args="$args --metadata_file $deseq_groups_file"
    echo "Expression-based prioritization enabled"
else
    args="$args --no_prioritize_by_expression"
    echo "Expression-based prioritization disabled"
fi

echo "Running module overview generation..."
# Check host script first (allows updates without rebuilding container)
if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/module_overview_interactive.py" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/module_overview_interactive.py"
elif [ -f "/app/scripts/module_overview_interactive.py" ]; then
    SCRIPT_PATH="/app/scripts/module_overview_interactive.py"
else
    echo "Error: module_overview_interactive.py not found"
    exit 1
fi
echo "Command: python3 $SCRIPT_PATH $args"

python3 $SCRIPT_PATH $args

echo "Module overview generation completed"
