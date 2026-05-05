#!/bin/bash -euo pipefail
# Stage all LemonPreprocessed_*.txt files (expression, complete, omics-specific)
    mkdir -p Preprocessing
    for f in LemonPreprocessed_complete.txt LemonPreprocessed_expression.txt LemonPreprocessed_metabolomics.txt; do
        if [ -f "$f" ]; then
            cp -L "$f" Preprocessing/ || true
        fi
    done

    mkdir -p heatmaps
    # Ensure canonical ModuleViewer_files directory exists and copy staged viewer files into it
    mkdir -p ModuleViewer_files
    for f in Metabolite_Gene_enrichment_results.csv Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt PPI_enrichment_results.csv PPI_interactions.mvf TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt evaluation_summary.txt metabolite_LemoniteKG_interactions.mvf sample_mapping.mvf; do
        # follow symlinks (-L) to ensure files referenced by symlinks are copied
        if [ -d "$f" ]; then
            cp -L -r "$f"/* ModuleViewer_files/ || true
        else
            cp -L -r "$f" ModuleViewer_files/ || true
        fi
    done
    # sample_mapping.mvf is expected to be present inside the staged viewer_files directory
    
    # Also stage preprocessing file(s) so the module viewer can find expression/omics data.
    # Need LemonPreprocessed_expression.txt (main heatmap), LemonPreprocessed_complete.txt (TF regulators),
    # and LemonPreprocessed_*.txt for each omics type (e.g. LemonPreprocessed_proteomics.txt, LemonPreprocessed_metabolomics.txt).
    # Try a few strategies and follow symlinks: 1) look under the staged input_dir, 2) search the workdir, 3) fallback to results paths.
    mkdir -p Preprocessing
    
    # Copy ALL LemonPreprocessed_*.txt files (expression, complete, and omics-specific)
    # 1) staged input_dir paths (if available)
    if [ -d "test_dataset/Preprocessing" ]; then
        for f in test_dataset/Preprocessing/LemonPreprocessed_*.txt; do
            [ -f "$f" ] && cp -L "$f" Preprocessing/ || true
        done
    fi
    for f in test_dataset/LemonPreprocessed_*.txt; do
        [ -f "$f" ] && cp -L "$f" Preprocessing/ || true
    done
    # 2) search the current workdir (shallow) for any LemonPreprocessed_*.txt
    # Use -maxdepth 2 to avoid following the input_dir symlink into work/ (causes a loop)
    find -L . -maxdepth 2 -type f -name 'LemonPreprocessed_*.txt' -print 2>/dev/null | while read -r src; do
        fname=$(basename "$src")
        [ ! -f "Preprocessing/$fname" ] && cp -L "$src" Preprocessing/ || true
    done || true
    # 3) try common results path(s) if present
    for p in ./results/*/LemonTree/Preprocessing ./results/*/Preprocessing ../results/*/LemonTree/Preprocessing; do
        for f in $p/LemonPreprocessed_*.txt; do
            if [ -f "$f" ]; then
                fname=$(basename "$f")
                [ ! -f "Preprocessing/$fname" ] && cp -L "$f" Preprocessing/ || true
            fi
        done
    done
    
    echo "Staged preprocessing files:"
    ls -la Preprocessing/LemonPreprocessed_*.txt 2>/dev/null || echo "  (none found)"
    
    # Build regulator_files argument dynamically based on regulator_types
    # Parse regulator_types (format: "Type:Prefix,Type:Prefix")
    IFS=',' read -ra REG_TYPES <<< "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"
    
    # Initialize regulator file variables
    declare -A regulator_files=()
    
    for reg_type in "${REG_TYPES[@]}"; do
        # Extract type and prefix (e.g., "TFs:lovering_TF_list.txt" -> type="TFs", prefix="lovering_TF_list.txt")
        IFS=':' read -r reg_type_name reg_prefix <<< "$reg_type"
        reg_type_name=$(echo $reg_type_name | xargs)  # trim whitespace
        reg_prefix=$(echo $reg_prefix | xargs)        # trim whitespace
        
        # Use the type name (Prefix) for the filename, NOT the data filename
        # Files are named like: TFs.selected_regs_list.txt, Metabolites.selected_regs_list.txt
        reg_file="ModuleViewer_files/${reg_type_name}.selected_regs_list.txt"
        if [ -f "$reg_file" ]; then
            regulator_files[$reg_type_name]="$reg_file"
            echo "Found ${reg_type_name} regulator file: $reg_file"
        else
            echo "Warning: ${reg_type_name} regulator file not found: $reg_file"
        fi
    done
    
    # Build regulator_files argument
    regulator_files_arg=""
    for reg_type in "${!regulator_files[@]}"; do
        reg_file="${regulator_files[$reg_type]}"
        if [ -n "$regulator_files_arg" ]; then
            regulator_files_arg="$regulator_files_arg,$reg_type:$reg_file"
        else
            regulator_files_arg="$reg_type:$reg_file"
        fi
    done
    
    # Run module viewer using the staged ModuleViewer_files (current working directory)
    # Check host script first (allows updates without rebuilding container)
    if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/module_viewer.py" ]; then
        SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/module_viewer.py"
    elif [ -f "/app/scripts/module_viewer.py" ]; then
        SCRIPT_PATH="/app/scripts/module_viewer.py"
    else
        echo "Error: module_viewer.py not found"
        exit 1
    fi
    PYTHONNOUSERSITE=1 python3 $SCRIPT_PATH \
        --input_dir . \
        --output_dir heatmaps \
        --regulator_files "$regulator_files_arg" \
        --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
        --dpi 300 \
        --show_regulator_scores \
        --annotation_types "diagnosis" \
        --name_mapping name_mapping.tsv \
        --modules $(cat specific_modules.txt | tr '
' ',' | sed 's/,$//')
    
    # Create summary
    echo "ModuleViewer Heatmaps Generated" > module_viewer_summary.txt
    echo "================================" >> module_viewer_summary.txt
    echo "Date: $(date)" >> module_viewer_summary.txt
    echo "Input directory: test_dataset" >> module_viewer_summary.txt
    echo "Viewer files: Metabolite_Gene_enrichment_results.csv Metabolites.selected_regs_list.txt Metabolites.selected_regulators_scores.txt PPI_enrichment_results.csv PPI_interactions.mvf TFs.selected_regs_list.txt TFs.selected_regulators_scores.txt clusters_list.txt evaluation_summary.txt metabolite_LemoniteKG_interactions.mvf sample_mapping.mvf" >> module_viewer_summary.txt
    echo "Regulator types: TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" >> module_viewer_summary.txt
    echo "" >> module_viewer_summary.txt
    echo "Generated files:" >> module_viewer_summary.txt
    ls -la heatmaps/ >> module_viewer_summary.txt
    echo "" >> module_viewer_summary.txt
    echo "Total heatmaps: $(ls heatmaps/*.png 2>/dev/null | wc -l)" >> module_viewer_summary.txt
