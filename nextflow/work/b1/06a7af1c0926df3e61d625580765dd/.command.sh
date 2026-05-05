#!/bin/bash -euo pipefail
# Create Networks directory and copy network files
mkdir -p Networks

# Copy all regulator target files to Networks directory (dynamically derived from --regulator_types)
for f in Metabolites2targets_top2.0pct_75_modules.txt TFs2targets_top2.0pct_75_modules.txt; do
    if [ -f "$f" ]; then
        cp "$f" Networks/
    fi
done

# Copy LemonNetwork files instead of Cytoscape files
for f in LemonNetwork_*.txt; do
    if [ -f "$f" ]; then
        cp "$f" Networks/ 2>/dev/null || true
    fi
done

# Copy preprocessed data files to working directory (avoid copying file onto itself)
if [ -f "LemonPreprocessed_expression.txt" ]; then
    # Check if source and destination are the same file to avoid copying onto itself
    if [ -f ./LemonPreprocessed_expression.txt ] && [ "LemonPreprocessed_expression.txt" -ef ./LemonPreprocessed_expression.txt ] 2>/dev/null; then
        echo "Preprocessed file already present at destination; skipping copy"
    else
        cp "LemonPreprocessed_expression.txt" ./LemonPreprocessed_expression.txt
    fi
fi

# Copy clusters file from viewer_files
if [ -d "Metabolite_Gene_enrichment_results.csv Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt PPI_enrichment_results.csv PPI_interactions.mvf TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt evaluation_summary.txt metabolite_LemoniteKG_interactions.mvf sample_mapping.mvf" ]; then
    # Look for clusters_list.txt in viewer_files directory
    find "Metabolite_Gene_enrichment_results.csv Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt PPI_enrichment_results.csv PPI_interactions.mvf TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt evaluation_summary.txt metabolite_LemoniteKG_interactions.mvf sample_mapping.mvf" -name "clusters_list.txt" -exec cp {} ./clusters_list.txt \; 2>/dev/null || true
fi

# Copy regulator list files from viewer_files
if [ -d "Metabolite_Gene_enrichment_results.csv Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt PPI_enrichment_results.csv PPI_interactions.mvf TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt evaluation_summary.txt metabolite_LemoniteKG_interactions.mvf sample_mapping.mvf" ]; then
    cp "Metabolite_Gene_enrichment_results.csv Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt PPI_enrichment_results.csv PPI_interactions.mvf TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt evaluation_summary.txt metabolite_LemoniteKG_interactions.mvf sample_mapping.mvf"/*_list.txt ./ 2>/dev/null || true
fi

# List what we have for debugging
echo "=== Files in working directory ==="
ls -la
echo "=== Files in Networks directory ==="
ls -la Networks/
echo "=================================="

# Run enrichment analysis
# Check host script first (allows updates without rebuilding container)
if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/enrichment_analysis.R" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/enrichment_analysis.R"
elif [ -f "/app/scripts/enrichment_analysis.R" ]; then
    SCRIPT_PATH="/app/scripts/enrichment_analysis.R"
else
    echo "Error: enrichment_analysis.R not found"
    exit 1
fi
Rscript $SCRIPT_PATH         --input_dir .         --output_dir .         --analysis_method "EnrichR"         --top_n_percent_regulators "2.0"         --coherence_threshold "0.5"         --n_threads 5         --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"         --organism "human"
