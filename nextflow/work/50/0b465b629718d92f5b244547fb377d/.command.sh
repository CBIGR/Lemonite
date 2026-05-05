#!/bin/bash -euo pipefail
# Resolve expression file: explicit param or auto-detect with uniqueness check
if [ -n "" ]; then
    EXPRESSION_FILE="null"
else
    MATCHES=$(find data/ -name "*host_tx_counts.tsv" -o -name "*expression*.tsv" -o -name "*counts*.tsv" -o -name "Counts.tsv")
    MATCH_COUNT=$(echo "$MATCHES" | grep -c . || true)
    if [ "$MATCH_COUNT" -eq 0 ]; then
        echo "Error: Expression file not found. Provide --expression_file or place a *counts*.tsv / Counts.tsv in data/"
        find data/ -type f -name "*.tsv" -o -name "*.txt"
        exit 1
    elif [ "$MATCH_COUNT" -gt 1 ]; then
        echo "Error: Multiple expression files matched ($MATCH_COUNT). Please specify --expression_file explicitly:"
        echo "$MATCHES"
        exit 1
    fi
    EXPRESSION_FILE=$(echo "$MATCHES" | head -1)
fi

# Resolve metadata file: explicit param or auto-detect with uniqueness check
if [ -n "" ]; then
    METADATA_FILE="null"
else
    MATCHES=$(find data/ -name "*metadata*.txt" -o -name "*Metadata*.txt")
    MATCH_COUNT=$(echo "$MATCHES" | grep -c . || true)
    if [ "$MATCH_COUNT" -eq 0 ]; then
        echo "Error: Metadata file not found. Provide --metadata_file or place a *metadata*.txt in data/"
        exit 1
    elif [ "$MATCH_COUNT" -gt 1 ]; then
        echo "Error: Multiple metadata files matched ($MATCH_COUNT). Please specify --metadata_file explicitly:"
        echo "$MATCHES"
        exit 1
    fi
    METADATA_FILE=$(echo "$MATCHES" | head -1)
fi

# Resolve prior network (optional)
PRIOR_NETWORK=$(find data/ -name "*CollecTRI*.txt" -o -name "*network*.txt" | head -1)
if [ -z "$PRIOR_NETWORK" ] && [ -f "/home/borisvdm/repo/LemonIte/nextflow/PKN/CollecTRI_network.txt" ]; then
    PRIOR_NETWORK="/home/borisvdm/repo/LemonIte/nextflow/PKN/CollecTRI_network.txt"
    echo "Using bundled network from projectDir: $PRIOR_NETWORK"
fi

echo "Found core files:"
echo "Expression: $EXPRESSION_FILE"
echo "Metadata: $METADATA_FILE"
echo "Prior network: $PRIOR_NETWORK"

# Copy name_map.csv if it exists
if [ -f "data/name_map.csv" ]; then
    mkdir -p "LemonTree/Preprocessing/"
    cp "data/name_map.csv" "LemonTree/Preprocessing/"
    echo "Copied name_map.csv to Preprocessing directory"
fi

# Run R preprocessing script with regulator_types parameter
# Check host script first (allows updates without rebuilding container)
if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/preprocessing_tfa_complete.R" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/preprocessing_tfa_complete.R"
elif [ -f "/app/scripts/preprocessing_tfa_complete.R" ]; then
    SCRIPT_PATH="/app/scripts/preprocessing_tfa_complete.R"
else
    echo "Error: preprocessing_tfa_complete.R not found"
    exit 1
fi
# Export gene annotation file path if provided; a local PKN backup will be available if online annotation fails.
export GENE_ANNOTATION_FILE=""
export GENE_ANNOTATION_BACKUP="/home/borisvdm/repo/LemonIte/nextflow/PKN/ensembl_mapping_jan2024.txt"

Rscript $SCRIPT_PATH \
    --expression "$EXPRESSION_FILE" \
    --metadata "$METADATA_FILE" \
    --output_dir . \
    --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
    --top_n_genes "2000"         --perform_TFA "true"         --use_omics_specific_scaling "true"         --DESeq_contrast1 "diagnosis"         --design_formula "~ diagnosis"         --metadata_columns "diagnosis"         --expression_col "count"         --sample_id_col "Sample_ID"         --organism "human"         $([ ! -z "$PRIOR_NETWORK" ] && echo "--prior_network $PRIOR_NETWORK")
