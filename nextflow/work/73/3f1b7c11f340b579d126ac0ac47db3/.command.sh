#!/bin/bash -euo pipefail
# Resolve expression file: explicit param or auto-detect with uniqueness check
if [ -n "/home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/data/proteomics_scaled.csv" ]; then
    EXPRESSION_FILE="/home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/data/proteomics_scaled.csv"
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
if [ -n "/home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/data/Metadata_AML_nextflow.csv" ]; then
    METADATA_FILE="/home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/data/Metadata_AML_nextflow.csv"
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

# Determine which preprocessing script to use
# "proteomics" mode: uses Preprocessing_TFA_Proteomics.R (pre-scaled, Pareto/z-score, TFA optional)
# "rna" mode (default): uses Preprocessing_TFA_RNA.R (DESeq2 normalisation, TFA enabled)
if [ "proteomics" = "proteomics" ] && [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/Preprocessing_TFA_Proteomics.R" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/Preprocessing_TFA_Proteomics.R"
    echo "Using Preprocessing_TFA_Proteomics.R (pre-scaled proteomics mode)"

    # Build optional file arguments
    OPTIONAL_ARGS=""

    # Check for ID mapping file
    if [ -f "data/Proteins_ID_mapping.tsv" ]; then
        OPTIONAL_ARGS="$OPTIONAL_ARGS --id_mapping data/Proteins_ID_mapping.tsv"
    fi

    # Parse regulator_types for hPTM and metabolomics files
    # Format: "hPTMs:Histone_proteomics.csv,Metabolites:metabolomics_combined.csv"
    REGULATOR_TYPES="hPTMs:Histone_proteomics.csv,Metabolites:metabolomics_combined.csv"
    if [[ "$REGULATOR_TYPES" == *"hPTM"* ]]; then
        HPTM_FILE=$(echo "$REGULATOR_TYPES" | grep -oP '(?<=hPTMs:)[^,]+')
        if [ -f "data/$HPTM_FILE" ]; then
            OPTIONAL_ARGS="$OPTIONAL_ARGS --hptm_file data/$HPTM_FILE"
        fi
    fi
    if [[ "$REGULATOR_TYPES" == *"Metabolites"* ]]; then
        METABO_FILE=$(echo "$REGULATOR_TYPES" | grep -oP '(?<=Metabolites:)[^,:]+')
        if [ -f "data/$METABO_FILE" ]; then
            OPTIONAL_ARGS="$OPTIONAL_ARGS --metabolomics_file data/$METABO_FILE"
        fi
        # Check for metabolomics labels file (explicit param takes priority, then data/ fallback)
        if [ -n "/home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/results/LemonTree/Proteomics/Preprocessing/Metabolites_name_map.txt" ] && [ -f "/home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/results/LemonTree/Proteomics/Preprocessing/Metabolites_name_map.txt" ]; then
            OPTIONAL_ARGS="$OPTIONAL_ARGS --metabolomics_labels '/home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/results/LemonTree/Proteomics/Preprocessing/Metabolites_name_map.txt'"
        elif [ -f "data/metabolomics_name_map.csv" ]; then
            OPTIONAL_ARGS="$OPTIONAL_ARGS --metabolomics_labels data/metabolomics_name_map.csv"
        fi
    fi

    Rscript $SCRIPT_PATH \
        --expression "$EXPRESSION_FILE" \
        --metadata "$METADATA_FILE" \
        --output_dir . \
        --top_n_genes "1000" \
        --perform_TFA "false" \
        --organism "human" \
        --sample_id_col "ID" \
        $OPTIONAL_ARGS
else
    # Default: RNA-seq preprocessing with DESeq2 and TFA
    if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/Preprocessing_TFA_RNA.R" ]; then
        SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/Preprocessing_TFA_RNA.R"
    elif [ -f "/app/scripts/Preprocessing_TFA_RNA.R" ]; then
        SCRIPT_PATH="/app/scripts/Preprocessing_TFA_RNA.R"
    else
        echo "Error: Neither preprocessing script found (Preprocessing_TFA_RNA.R)"
        exit 1
    fi
    echo "Using Preprocessing_TFA_RNA.R for RNA-seq count data"

    # Export gene annotation file path if provided
    export GENE_ANNOTATION_FILE=""
    export GENE_ANNOTATION_BACKUP="/home/borisvdm/repo/LemonIte/nextflow/PKN/ensembl_mapping_jan2024.txt"

    Rscript $SCRIPT_PATH \
        --expression "$EXPRESSION_FILE" \
        --metadata "$METADATA_FILE" \
        --output_dir . \
        --regulator_types "hPTMs:Histone_proteomics.csv,Metabolites:metabolomics_combined.csv" \
        --top_n_genes "1000"             --perform_TFA "false"             --use_omics_specific_scaling "true"             --DESeq_contrast1 "Cell_line"             --design_formula "~ Cell_line"             --metadata_columns "Cell_line,FAB,MITO_subtype,Risk_classification"             --expression_col "count"             --sample_id_col "ID"             --organism "human"             $([ ! -z "$PRIOR_NETWORK" ] && echo "--prior_network $PRIOR_NETWORK")
fi
