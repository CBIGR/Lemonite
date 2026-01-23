process MODULE_OVERVIEW_INTERACTIVE {
    tag "Interactive module overview"
    publishDir "${(params.output_dir ?: (params.input_dir ? params.input_dir + '/results' : './results'))}/${(params.computed_run_id ?: (params.run_id ?: 'run_auto'))}/LemonTree", mode: 'copy'

    input:
    path viewer_files
    path filtered_modules
    path coherence_scores
    path lemontree_outputs
    path expression_file
    path deseq_groups_file
    val run_id
    val regulator_types

    output:
    path "Module_Overview/*", emit: overview_results
    path "Module_Overview/Module_Expression_Heatmap.png", emit: expression_heatmap, optional: true

    script:
    def top_n_percent_regulators = params.top_n_percent_regulators ?: 2.0
    def coherence_threshold = params.coherence_threshold ?: 0.6
    def group_column = params.deseq_contrast1 ?: 'diagnosis'
    def n_clusters = params.overview_n_clusters ?: 5
    // Normalize enrichment method to proper case (enrichr -> EnrichR, gsea -> GSEA)
    def enrichment_input = params.enrichment_method ?: 'auto'
    def enrichment_method = enrichment_input.toLowerCase() == 'enrichr' ? 'EnrichR' : 
                           enrichment_input.toLowerCase() == 'gsea' ? 'GSEA' : 
                           enrichment_input
    def use_megago = params.use_megago
    def prioritize_by_expression = params.prioritize_by_expression
    
    """
    # Create working directory
    workdir=\$(pwd)
    
    # Find required files
    cluster_file=\$(find . -name "clusters_list.txt" | head -1)
    coherence_file="${coherence_scores}"
    
    # Use the provided expression and DESeq files
    expression_file="${expression_file}"
    deseq_groups_file="${deseq_groups_file}"
    
    # Also look for the MVF file that should be created by PKN evaluation
    mvf_file=\$(find . -name "metabolite_LemoniteKG_interactions.mvf" | head -1)
    
    echo "=== Module Overview Generation ==="
    echo "Found files:"
    echo "  Clusters: \$cluster_file"
    echo "  Coherence scores: \$coherence_file"
    echo "Expression: \$expression_file (provided as input)"
    echo "  DESeq groups: \$deseq_groups_file (provided as input)"
    echo "  MVF file: \$mvf_file"
    echo "Regulator types: ${regulator_types}"
    echo "Parameters:"
    echo "  Enrichment method: ${enrichment_method}"
    echo "  N clusters: ${n_clusters}"
    echo "  Coherence threshold: ${coherence_threshold}"
    echo "  Use megago: ${use_megago}"
    # Debug: list enrichment dir if present
    if [ -d "Enrichment" ]; then
        echo "=== Enrichment directory contents ==="
        ls -la Enrichment || true
        echo "====================================="
    fi
    
    # Build regulator_files argument dynamically based on regulator_types
    REGULATOR_FILES=""
    REGULATOR_TYPES="${regulator_types}"
    
    # Parse regulator types and find corresponding files
    IFS=',' read -ra TYPE_ARRAY <<< "\$REGULATOR_TYPES"
    for type_config in "\${TYPE_ARRAY[@]}"; do
        # Extract type name (before colon) - this is the Prefix like "TFs", "Metabolites"
        TYPE_NAME=\$(echo "\$type_config" | cut -d':' -f1)
        # Extract data filename (after colon) - not used for file lookup
        FILE_PREFIX=\$(echo "\$type_config" | cut -d':' -f2)
        
        # Find the corresponding regulator file using TYPE_NAME (Prefix), NOT data filename
        # Files are named like: TFs.selected_regs_list.txt, Metabolites.selected_regs_list.txt
        REG_FILE="\${TYPE_NAME}.selected_regs_list.txt"
        if [ -f "\$REG_FILE" ]; then
            if [ -n "\$REGULATOR_FILES" ]; then
                REGULATOR_FILES="\$REGULATOR_FILES,\${TYPE_NAME}:\${REG_FILE}"
            else
                REGULATOR_FILES="\${TYPE_NAME}:\${REG_FILE}"
            fi
            echo "Found \${TYPE_NAME} regulator file: \$REG_FILE"
        else
            echo "Warning: Regulator file \$REG_FILE not found"
        fi
    done
    
    # Construct command arguments
    args="--input_dir . --output_dir ."
    args="\$args --enrichment_method ${enrichment_method}"
    args="\$args --coherence_threshold ${coherence_threshold}"
    args="\$args --group_column ${group_column}"
    if [ -n "\$REGULATOR_FILES" ]; then
        args="\$args --regulator_files \$REGULATOR_FILES"
    fi
    
    # Add megago clustering if enabled
    if [ "${use_megago}" = "true" ]; then
        args="\$args --n_clusters ${n_clusters}"
        echo "Megago functional clustering enabled with ${n_clusters} clusters"
    else
        args="\$args --n_clusters 1"
        echo "Megago clustering disabled - using simple module grouping"
    fi
    
    # Expression prioritization is enabled by default but can be disabled
    if [ "${prioritize_by_expression}" = "true" ]; then
        args="\$args --prioritize_by_expression"
        args="\$args --expression_file \$expression_file"
        args="\$args --metadata_file \$deseq_groups_file"
        echo "Expression-based prioritization enabled"
    else
        args="\$args --no_prioritize_by_expression"
        echo "Expression-based prioritization disabled"
    fi
    
    echo "Running module overview generation..."
    # Check host script first (allows updates without rebuilding container)
    if [ -f "${projectDir}/scripts/module_overview_interactive.py" ]; then
        SCRIPT_PATH="${projectDir}/scripts/module_overview_interactive.py"
    elif [ -f "/app/scripts/module_overview_interactive.py" ]; then
        SCRIPT_PATH="/app/scripts/module_overview_interactive.py"
    else
        echo "Error: module_overview_interactive.py not found"
        exit 1
    fi
    echo "Command: python3 \$SCRIPT_PATH \$args"
    
    python3 \$SCRIPT_PATH \$args
    
    echo "Module overview generation completed"
    """

    stub:
    """
    mkdir -p Module_Overview
    touch Module_Overview/Module_Overview_Interactive.csv
    touch Module_Overview/module_functional_clusters.csv
    touch Module_Overview/interactive_module_network.html
    touch Module_Overview/cluster_characteristics_heatmap.html
    touch Module_Overview/cluster_statistics.csv
    touch Module_Overview/module_expression_analysis.csv
    """
}
