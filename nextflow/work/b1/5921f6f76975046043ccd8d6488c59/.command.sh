#!/bin/bash -euo pipefail
# Create directory structure for report generation
mkdir -p LemonTree/Networks LemonTree/ModuleViewer_files LemonTree/Enrichment LemonTree/PKN_Evaluation LemonTree/Preprocessing

echo "=== Files staged by Nextflow ==="
ls -la
echo "================================="

# Copy viewer files
for f in Metabolite_Gene_enrichment_results.csv Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt PPI_enrichment_results.csv PPI_interactions.mvf TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt evaluation_summary.txt metabolite_LemoniteKG_interactions.mvf sample_mapping.mvf; do
    if [ -d "$f" ]; then
        cp -rL "$f"/* LemonTree/ModuleViewer_files/ 2>/dev/null || true
    else
        cp -L "$f" LemonTree/ModuleViewer_files/ 2>/dev/null || true
    fi
done

# Copy network files
for f in Cytoscape_attributes_75_modules_filtered.txt Cytoscape_network_75_modules_top2.0pct_filtered.txt; do
    if [ -d "$f" ]; then
        cp -rL "$f"/* LemonTree/Networks/ 2>/dev/null || true
    else
        cp -L "$f" LemonTree/Networks/ 2>/dev/null || true
    fi
done

# Copy filtered modules
if [ -f "specific_modules.txt" ]; then
    cp -L "specific_modules.txt" LemonTree/Networks/specific_modules.txt 2>/dev/null || true
fi

# Copy coherence scores file
if [ -f "Module_coherence_scores.txt" ]; then
    cp -L "Module_coherence_scores.txt" LemonTree/Networks/Module_coherence_scores.txt 2>/dev/null || true
fi

# Copy enrichment results - skip the large PNG directories to avoid memory issues
if [ -d "Enrichment" ]; then
    # Only copy CSV files, not PNG images
    find "Enrichment" -maxdepth 3 -name "*.csv" -exec cp --parents {} LemonTree/Enrichment/ \; 2>/dev/null || true
elif [ -f "Enrichment" ]; then
    cp -L "Enrichment" LemonTree/Enrichment/ 2>/dev/null || true
fi

# Copy PKN results
for f in ; do
    if [ -d "$f" ]; then
        cp -rL "$f"/* LemonTree/PKN_Evaluation/ 2>/dev/null || true
    elif [ -f "$f" ]; then
        cp -L "$f" LemonTree/PKN_Evaluation/ 2>/dev/null || true
    fi
done

# Copy preprocessing results
for f in DESeq_groups.txt LemonPreprocessed_complete.txt LemonPreprocessed_expression.txt LemonPreprocessed_metabolomics.txt Normalized_metabolites.pdf PCA_metabolites.pdf PCA_transcriptomics.pdf expression_variance_histogram.pdf metabolites.txt name_map.csv name_mapping.tsv tfs.txt; do
    if [ -d "$f" ]; then
        cp -rL "$f"/* LemonTree/Preprocessing/ 2>/dev/null || true
    elif [ -f "$f" ]; then
        cp -L "$f" LemonTree/Preprocessing/ 2>/dev/null || true
    fi
done

# Copy clustering results (Lemon_out directory with tight_clusters and cluster dirs)
mkdir -p LemonTree/Lemon_out
for f in cluster_3 cluster_2 cluster_1 cluster_4 cluster_5; do
    if [ -d "$f" ]; then
        # If it's a directory, copy its contents
        cp -rL "$f"/* LemonTree/Lemon_out/ 2>/dev/null || true
    elif [ -f "$f" ]; then
        # If it's a file, copy it directly
        cp -L "$f" LemonTree/Lemon_out/ 2>/dev/null || true
    fi
done

# Copy Module_Overview results (Module_Overview.csv, module_expression_analysis.csv, etc.)
mkdir -p LemonTree/Module_Overview
for f in Module_Overview.csv Module_Overview_Comprehensive.html interactive_module_network.html interactive_module_network_movable.html module_expression_analysis.csv module_network_edges.txt module_network_node_attributes.txt regulator_rankings.html top_30; do
    if [ -d "$f" ]; then
        # If it's a directory, copy its contents
        cp -rL "$f"/* LemonTree/Module_Overview/ 2>/dev/null || true
    elif [ -f "$f" ]; then
        # If it's a file, copy it directly
        cp -L "$f" LemonTree/Module_Overview/ 2>/dev/null || true
    fi
done

# Copy parameters log
if [ -f "pipeline_parameters_log.txt" ]; then
    cp -L "pipeline_parameters_log.txt" ./pipeline_parameters_log.txt 2>/dev/null || true
fi

# Debug: list what we have
echo "=== LemonTree directory structure ==="
find LemonTree -type f | head -50 || true
echo "======================================"

# Run summary report generator
if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/generate_summary_report.py" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/generate_summary_report.py"
elif [ -f "/app/scripts/generate_summary_report.py" ]; then
    SCRIPT_PATH="/app/scripts/generate_summary_report.py"
else
    echo "Error: generate_summary_report.py not found"
    exit 1
fi

# Build name mapping argument — find it from staged preprocessing results
NAME_MAPPING_ARG=""
NAME_MAPPING_FILE=$(find ./LemonTree/Preprocessing -name "name_mapping.tsv" 2>/dev/null | head -1)
if [ -n "${NAME_MAPPING_FILE}" ] && [ -f "${NAME_MAPPING_FILE}" ]; then
    NAME_MAPPING_ARG="--name_mapping ${NAME_MAPPING_FILE}"
fi

# Pass current directory as output_dir with run_id="." so script looks for LemonTree directly
# The LemonTree directory is already copied to the current working directory by this point
python3 $SCRIPT_PATH \
    --input_dir "/home/borisvdm/repo/LemonIte/nextflow/test_dataset" \
    --output_dir . \
    --run_id . \
    --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
    --parameters_file ./pipeline_parameters_log.txt \
    --organism "human" \
    $NAME_MAPPING_ARG

# Move the report to the expected output location
if [ -f "2000HVG_coherence0.5_top2.0pct_clusters5/Lemonite_Summary_Report.html" ]; then
    mv "2000HVG_coherence0.5_top2.0pct_clusters5/Lemonite_Summary_Report.html" ./Lemonite_Summary_Report.html
fi

echo "Summary report generation completed"
