#!/bin/bash -euo pipefail
# Create directory structure for network generation
mkdir -p Lemon_out Preprocessing

echo "=== Files staged by Nextflow ==="
ls -la
echo "================================="

# Copy all lemontree_outputs to Lemon_out directory
# These are staged directly in the work directory by Nextflow
# Use -L to follow symlinks and copy actual file contents
for file in Lemon_results Metabolites.allreg.txt Metabolites.randomreg.txt Metabolites.topreg.txt Metabolites.xml.gz clusterfile hPTMs.allreg.txt hPTMs.randomreg.txt hPTMs.topreg.txt hPTMs.xml.gz tight_clusters.txt; do
    if [ -L "$file" ]; then
        # If it's a symlink, copy the target
        cp -rL "$file" ./Lemon_out/
    elif [ -f "$file" ]; then
        # If it's a regular file, copy it directly
        cp "$file" ./Lemon_out/
    elif [ -d "$file" ]; then
        # If it's a directory, copy it recursively (dereference symlinks for HPC compatibility)
        cp -rL "$file" ./Lemon_out/
    fi
done

# Copy all preprocessing_files to Preprocessing directory
for file in DESeq_groups.txt LemonPreprocessed_complete.txt LemonPreprocessed_expression.txt LemonPreprocessed_hPTM.txt LemonPreprocessed_metabolomics.txt LemonPreprocessed_proteomics.txt PCA_hptm.pdf PCA_metabolomics.pdf PCA_proteomics.pdf hptms.txt metabolites.txt name_mapping.tsv; do
    if [ -L "$file" ]; then
        cp -L "$file" ./Preprocessing/
    else
        cp "$file" ./Preprocessing/
    fi
done

# List available files for debugging
echo "=== Available files before network generation ==="
echo "=== Lemon_out contents ==="
ls -la Lemon_out/
find Lemon_out/ -name "*.allreg.txt" -o -name "*.topreg.txt"
echo "=== Preprocessing contents ==="
ls -la Preprocessing/
echo "================================================="

# Validate key input files before network generation
if [ ! -d "Lemon_out" ] || [ -z "$(ls -A Lemon_out/)" ]; then
    echo "Error: Lemon_out directory is missing or empty"
    exit 1
fi
if [ ! -d "Preprocessing" ] || [ ! -f "Preprocessing/LemonPreprocessed_expression.txt" ]; then
    echo "Error: Preprocessing directory or expression file missing"
    exit 1
fi

# Run network generation script
echo "Running network generation..."
# Check host script first (allows updates without rebuilding container)
if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/lemontree_to_network.py" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/lemontree_to_network.py"
elif [ -f "/app/scripts/lemontree_to_network.py" ]; then
    SCRIPT_PATH="/app/scripts/lemontree_to_network.py"
else
    echo "Error: lemontree_to_network.py not found"
    exit 1
fi
python3 $SCRIPT_PATH         --input_dir .         --output_dir .         --run_id minW0p10         --coherence_threshold "0.5"         --regulator_selection_method "fold_per_module"         --top_n_percent_regulators "2.0"         --regulator_fold_cutoff "2.0"         --regulator_types "hPTMs:Histone_proteomics.csv,Metabolites:metabolomics_combined.csv"
