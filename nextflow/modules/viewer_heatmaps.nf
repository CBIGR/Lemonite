process MODULE_VIEWER_HEATMAPS {
    tag "Module viewer heatmaps"
    publishDir "${(params.output_dir ?: (params.input_dir ? params.input_dir + '/results' : './results'))}/${(params.computed_run_id ?: (params.run_id ?: 'run_auto'))}/LemonTree/module_heatmaps", mode: 'copy'
    
    input:
    path viewer_files
    path filtered_modules
    path input_dir
    path preprocessed_file
    val run_id
    val regulator_types
    
    output:
    path "heatmaps/*", emit: heatmaps
    path "module_viewer_summary.txt", emit: summary
    
    script:
    def modules_arg = filtered_modules.name != 'NO_FILE' ? "--modules \$(cat ${filtered_modules} | tr '\n' ',')" : ""
    def top_n_percent_regulators = params.top_n_percent_regulators ?: 2.0
    
    """
    # Create output directory
    mkdir -p heatmaps
    # Ensure canonical ModuleViewer_files directory exists and copy staged viewer files into it
    mkdir -p ModuleViewer_files
    for f in ${viewer_files}; do
        # follow symlinks (-L) to ensure files referenced by symlinks are copied
        if [ -d "\$f" ]; then
            cp -L -r "\$f"/* ModuleViewer_files/ || true
        else
            cp -L -r "\$f" ModuleViewer_files/ || true
        fi
    done
    # sample_mapping.mvf is expected to be present inside the staged viewer_files directory
    
    # Also stage preprocessing file(s) so the module viewer can find expression data.
    # Need both LemonPreprocessed_expression.txt (for main heatmap) and LemonPreprocessed_complete.txt (for regulator blocks)
    # Try a few strategies and follow symlinks: 1) look under the staged input_dir, 2) search the workdir (following symlinks), 3) fallback to any results/*/LemonTree/Preprocessing path.
    mkdir -p Preprocessing
    
    # Copy LemonPreprocessed_expression.txt
    # 1) staged input_dir paths (if available)
    if [ -d "${input_dir}/Preprocessing" ] && [ -f "${input_dir}/Preprocessing/LemonPreprocessed_expression.txt" ]; then
        cp -L "${input_dir}/Preprocessing/LemonPreprocessed_expression.txt" Preprocessing/ || true
    fi
    if [ -f "${input_dir}/LemonPreprocessed_expression.txt" ]; then
        cp -L "${input_dir}/LemonPreprocessed_expression.txt" Preprocessing/ || true
    fi
    # 2) search the current workdir and follow symlinks
    PREPROC_EXPR=\$(find -L . -type f -name 'LemonPreprocessed_expression.txt' -print -quit 2>/dev/null || true)
    if [ -n "\$PREPROC_EXPR" ]; then
        cp -L "\$PREPROC_EXPR" Preprocessing/ || true
    fi
    # 3) try common results path(s) if present
    for p in ./results/*/LemonTree/Preprocessing ./results/*/Preprocessing ../results/*/LemonTree/Preprocessing; do
        if ls \$p/LemonPreprocessed_expression.txt 1> /dev/null 2>&1; then
            cp -L \$p/LemonPreprocessed_expression.txt Preprocessing/ || true
            break
        fi
    done
    
    # Copy LemonPreprocessed_complete.txt
    # 1) staged input_dir paths (if available)
    if [ -d "${input_dir}/Preprocessing" ] && [ -f "${input_dir}/Preprocessing/LemonPreprocessed_complete.txt" ]; then
        cp -L "${input_dir}/Preprocessing/LemonPreprocessed_complete.txt" Preprocessing/ || true
    fi
    if [ -f "${input_dir}/LemonPreprocessed_complete.txt" ]; then
        cp -L "${input_dir}/LemonPreprocessed_complete.txt" Preprocessing/ || true
    fi
    # 1b) copy the explicitly passed preprocessed_file (deterministic channel)
    if [ -n "${preprocessed_file}" ] && [ -f "${preprocessed_file}" ]; then
        cp -L "${preprocessed_file}" Preprocessing/ || true
    fi
    # 2) search the current workdir and follow symlinks to find any LemonPreprocessed_complete.txt
    PREPROC_SRC=\$(find -L . -type f -name 'LemonPreprocessed_complete.txt' -print -quit 2>/dev/null || true)
    if [ -n "\$PREPROC_SRC" ]; then
        cp -L "\$PREPROC_SRC" Preprocessing/ || true
    fi
    # 3) try common results path(s) if present
    for p in ./results/*/LemonTree/Preprocessing ./results/*/Preprocessing ../results/*/LemonTree/Preprocessing; do
        if ls \$p/LemonPreprocessed_complete.txt 1> /dev/null 2>&1; then
            cp -L \$p/LemonPreprocessed_complete.txt Preprocessing/ || true
            break
        fi
    done
    
    # Build regulator_files argument dynamically based on regulator_types
    # Parse regulator_types (format: "Type:Prefix,Type:Prefix")
    IFS=',' read -ra REG_TYPES <<< "${regulator_types}"
    
    # Initialize regulator file variables
    declare -A regulator_files=()
    
    for reg_type in "\${REG_TYPES[@]}"; do
        # Extract type and prefix (e.g., "TFs:lovering_TF_list.txt" -> type="TFs", prefix="lovering_TF_list.txt")
        IFS=':' read -r reg_type_name reg_prefix <<< "\$reg_type"
        reg_type_name=\$(echo \$reg_type_name | xargs)  # trim whitespace
        reg_prefix=\$(echo \$reg_prefix | xargs)        # trim whitespace
        
        # Use the type name (Prefix) for the filename, NOT the data filename
        # Files are named like: TFs.selected_regs_list.txt, Metabolites.selected_regs_list.txt
        reg_file="ModuleViewer_files/\${reg_type_name}.selected_regs_list.txt"
        if [ -f "\$reg_file" ]; then
            regulator_files[\$reg_type_name]="\$reg_file"
            echo "Found \${reg_type_name} regulator file: \$reg_file"
        else
            echo "Warning: \${reg_type_name} regulator file not found: \$reg_file"
        fi
    done
    
    # Build regulator_files argument
    regulator_files_arg=""
    for reg_type in "\${!regulator_files[@]}"; do
        reg_file="\${regulator_files[\$reg_type]}"
        if [ -n "\$regulator_files_arg" ]; then
            regulator_files_arg="\$regulator_files_arg,\$reg_type:\$reg_file"
        else
            regulator_files_arg="\$reg_type:\$reg_file"
        fi
    done
    
    # Run module viewer using the staged ModuleViewer_files (current working directory)
    # Check host script first (allows updates without rebuilding container)
    if [ -f "${projectDir}/scripts/module_viewer.py" ]; then
        SCRIPT_PATH="${projectDir}/scripts/module_viewer.py"
    elif [ -f "/app/scripts/module_viewer.py" ]; then
        SCRIPT_PATH="/app/scripts/module_viewer.py"
    else
        echo "Error: module_viewer.py not found"
        exit 1
    fi
    python3 \$SCRIPT_PATH \\
        --input_dir . \\
        --output_dir heatmaps \\
        --regulator_files "\$regulator_files_arg" \\
        --dpi 300 \\
        --show_regulator_scores \\
        --annotation_types "${params.heatmap_metadata_cols ?: params.metadata_columns ?: 'diagnosis'}" \\
        ${modules_arg}
    
    # Create summary
    echo "ModuleViewer Heatmaps Generated" > module_viewer_summary.txt
    echo "================================" >> module_viewer_summary.txt
    echo "Date: \$(date)" >> module_viewer_summary.txt
    echo "Input directory: ${input_dir}" >> module_viewer_summary.txt
    echo "Viewer files: ${viewer_files}" >> module_viewer_summary.txt
    echo "Regulator types: ${regulator_types}" >> module_viewer_summary.txt
    echo "" >> module_viewer_summary.txt
    echo "Generated files:" >> module_viewer_summary.txt
    ls -la heatmaps/ >> module_viewer_summary.txt
    echo "" >> module_viewer_summary.txt
    echo "Total heatmaps: \$(ls heatmaps/*.png 2>/dev/null | wc -l)" >> module_viewer_summary.txt
    """
}
