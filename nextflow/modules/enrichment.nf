process ENRICHMENT_ANALYSIS {
    tag "Pathway enrichment analysis"
    publishDir "${(params.output_dir ?: (params.input_dir ? params.input_dir + '/results' : './results'))}/${(params.computed_run_id ?: (params.run_id ?: 'run_auto'))}/LemonTree", mode: 'copy'

    input:
    path filtered_modules
    path metabolites_targets
    path tf_targets
    path lipids_targets
    path network_files
    path viewer_files
    path preprocessed_data
    val run_id

    output:
    // Publish the entire Enrichment directory so downstream processes have a ./Enrichment folder
    path "Enrichment", emit: enrichment_results

    script:
    def analysis_method = params.enrichment_method ?: "EnrichR"
    if (analysis_method == "auto") {
        analysis_method = "both"
    }
    """
    # Create Networks directory and copy network files
    mkdir -p Networks
    
    # Copy target files to Networks directory
    for f in ${metabolites_targets} ${tf_targets} ${lipids_targets}; do
        if [ -f "\$f" ]; then
            cp "\$f" Networks/
        fi
    done
    
    # Copy LemonNetwork files instead of Cytoscape files
    for f in LemonNetwork_*.txt; do
        if [ -f "\$f" ]; then
            cp "\$f" Networks/ 2>/dev/null || true
        fi
    done
    
    # Copy preprocessed data files to working directory (avoid copying file onto itself)
    if [ -f "${preprocessed_data}" ]; then
        # Check if source and destination are the same file to avoid copying onto itself
        if [ -f ./LemonPreprocessed_expression.txt ] && [ "${preprocessed_data}" -ef ./LemonPreprocessed_expression.txt ] 2>/dev/null; then
            echo "Preprocessed file already present at destination; skipping copy"
        else
            cp "${preprocessed_data}" ./LemonPreprocessed_expression.txt
        fi
    fi
    
    # Copy clusters file from viewer_files
    if [ -d "${viewer_files}" ]; then
        # Look for clusters_list.txt in viewer_files directory
        find "${viewer_files}" -name "clusters_list.txt" -exec cp {} ./clusters_list.txt \\; 2>/dev/null || true
    fi
    
    # Copy regulator list files from viewer_files
    if [ -d "${viewer_files}" ]; then
        cp "${viewer_files}"/*_list.txt ./ 2>/dev/null || true
    fi
    
    # List what we have for debugging
    echo "=== Files in working directory ==="
    ls -la
    echo "=== Files in Networks directory ==="
    ls -la Networks/
    echo "=================================="
    
    # Run enrichment analysis
    # Check host script first (allows updates without rebuilding container)
    if [ -f "${projectDir}/scripts/enrichment_analysis.R" ]; then
        SCRIPT_PATH="${projectDir}/scripts/enrichment_analysis.R"
    elif [ -f "/app/scripts/enrichment_analysis.R" ]; then
        SCRIPT_PATH="/app/scripts/enrichment_analysis.R"
    else
        echo "Error: enrichment_analysis.R not found"
        exit 1
    fi
    Rscript \$SCRIPT_PATH \
        --input_dir . \
        --output_dir . \
        --analysis_method ${analysis_method} \
        --top_n_percent_regulators ${params.top_n_percent_regulators} \
        --coherence_threshold ${params.coherence_threshold} \
        --n_threads ${task.cpus} \
        --regulator_types "${params.regulator_types}" \
        --organism ${params.organism}
    """

    stub:
    """
    mkdir -p Enrichment/GSEA Enrichment/EnrichR Networks
    touch Enrichment/enrichment_summary.csv
    touch Networks/targets.txt
    """
}
