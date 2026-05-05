#!/bin/bash -euo pipefail
# Create directory structure
mkdir -p Lemon_out ModuleViewer_files Preprocessing Networks/subnetworks

echo "=== Files staged by Nextflow ==="
ls -la
echo "================================="

# Copy all viewer files to ModuleViewer_files
for file in Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt sample_mapping.mvf; do
    cp "$file" ./ModuleViewer_files/
done

# Copy all preprocessing files to Preprocessing
for file in DESeq_groups.txt LemonPreprocessed_complete.txt LemonPreprocessed_expression.txt LemonPreprocessed_metabolomics.txt Normalized_metabolites.pdf PCA_metabolites.pdf PCA_transcriptomics.pdf expression_variance_histogram.pdf metabolites.txt name_map.csv name_mapping.tsv tfs.txt; do
    cp "$file" ./Preprocessing/
done

echo "=== ModuleViewer_files contents ==="
ls -la ModuleViewer_files/
echo "=== Preprocessing contents ==="
ls -la Preprocessing/
echo "================================="

# Find all regulator files dynamically based on regulator_types parameter
REGULATOR_FILES=""
REGULATOR_TYPES="TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"

# Parse regulator types and find corresponding files
IFS=',' read -ra TYPE_ARRAY <<< "$REGULATOR_TYPES"
for type_config in "${TYPE_ARRAY[@]}"; do
    # Extract prefix (before colon)
    PREFIX=$(echo "$type_config" | cut -d':' -f1)

    # Find the corresponding regulator file using lowercase prefix
    PREFIX_LOWER=$(echo "$PREFIX" | tr '[:upper:]' '[:lower:]')
    REG_FILE=$(find ModuleViewer_files -name "${PREFIX}.*" -type f | head -1)
    if [ $(find ModuleViewer_files -name "${PREFIX}.*" -type f | wc -l) -gt 1 ]; then
        echo "Warning: Multiple files found for ${PREFIX}, using: ${REG_FILE}"
    fi

    if [ -n "$REG_FILE" ]; then
        if [ -z "$REGULATOR_FILES" ]; then
            REGULATOR_FILES="${PREFIX}:${REG_FILE}"
        else
            REGULATOR_FILES="${REGULATOR_FILES},${PREFIX}:${REG_FILE}"
        fi
        echo "Found ${PREFIX} regulators: ${REG_FILE}"
    else
        echo "WARNING: ${PREFIX} regulators file not found"
    fi
done

# Check if we have any regulator files
if [ -z "$REGULATOR_FILES" ]; then
    echo "WARNING: No regulator files found, skipping subnetwork graphs"
    exit 0
fi

# Find clusters file
CLUSTERS_FILE=$(find . -name "clusters_list.txt" -o -name "specific_modules.txt" | head -1)
if [ -z "$CLUSTERS_FILE" ] || [ ! -f "$CLUSTERS_FILE" ]; then
    echo "WARNING: Clusters file not found, skipping subnetwork graphs"
    exit 0
fi

# Copy clusters file to expected location
mkdir -p ./Lemon_out
cp "$CLUSTERS_FILE" ./Lemon_out/clusters_list.txt

# Run subnetwork graph generation
echo "Creating subnetwork visualization graphs..."
# Check host script first (allows updates without rebuilding container)
if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/create_subnetwork_graphs.py" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/create_subnetwork_graphs.py"
elif [ -f "/app/scripts/create_subnetwork_graphs.py" ]; then
    SCRIPT_PATH="/app/scripts/create_subnetwork_graphs.py"
else
    echo "Error: create_subnetwork_graphs.py not found"
    exit 1
fi

# Check if metabolite mapping exists and build argument conditionally
METABOLITE_ARG=""
if [ -f ./Preprocessing/name_map.csv ]; then
    METABOLITE_ARG="--metabolite_mapping ./Preprocessing/name_map.csv"
    echo "Using metabolite mapping file: ./Preprocessing/name_map.csv"
else
    echo "No metabolite mapping file found - proceeding without metabolite name mapping"
fi

# Check if name mapping exists for original name restoration
NAME_MAPPING_ARG=""
NAME_MAPPING_FILE=$(find ./Preprocessing -name "name_mapping.tsv" | head -1)
if [ -n "${NAME_MAPPING_FILE}" ] && [ -f "${NAME_MAPPING_FILE}" ]; then
    NAME_MAPPING_ARG="--name_mapping ${NAME_MAPPING_FILE}"
    echo "Using name mapping file for original name restoration: ${NAME_MAPPING_FILE}"
fi

python3 $SCRIPT_PATH         --regulator_files "$REGULATOR_FILES"         --clusters ./Lemon_out/clusters_list.txt         --pkn "/home/borisvdm/repo/LemonIte/nextflow/PKN/Lemonite_PKN.tsv"         $METABOLITE_ARG         $NAME_MAPPING_ARG         --output_dir ./Networks/subnetworks

echo "Subnetwork graphs created successfully"
