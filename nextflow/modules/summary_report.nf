process SUMMARY_REPORT {
    tag "Generating summary report"
    publishDir "${(params.output_dir ?: (params.input_dir ? params.input_dir + '/results' : './results'))}/${(params.computed_run_id ?: (params.run_id ?: 'run_auto'))}", mode: 'copy'
    
    memory '8 GB'
    cpus 2
    time '1h'

    input:
    path viewer_files
    path filtered_modules
    path network_files
    path coherence_scores
    path enrichment_results
    path pkn_results
    path preprocessing_results
    path clustering_results
    path overview_results
    path parameters_log
    val run_id
    val regulator_types

    output:
    path "Lemonite_Summary_Report.html", emit: summary_report
    path "LemonTree/Module_Overview/*", emit: module_overview_files, optional: true

    script:
    """
    # Create directory structure for report generation
    mkdir -p LemonTree/Networks LemonTree/ModuleViewer_files LemonTree/Enrichment LemonTree/PKN_Evaluation LemonTree/Preprocessing
    
    echo "=== Files staged by Nextflow ==="
    ls -la
    echo "================================="
    
    # Copy viewer files
    for f in ${viewer_files}; do
        if [ -d "\$f" ]; then
            cp -rL "\$f"/* LemonTree/ModuleViewer_files/ 2>/dev/null || true
        else
            cp -L "\$f" LemonTree/ModuleViewer_files/ 2>/dev/null || true
        fi
    done
    
    # Copy network files
    for f in ${network_files}; do
        if [ -d "\$f" ]; then
            cp -rL "\$f"/* LemonTree/Networks/ 2>/dev/null || true
        else
            cp -L "\$f" LemonTree/Networks/ 2>/dev/null || true
        fi
    done
    
    # Copy filtered modules
    if [ -f "${filtered_modules}" ]; then
        cp -L "${filtered_modules}" LemonTree/Networks/specific_modules.txt 2>/dev/null || true
    fi
    
    # Copy coherence scores file
    if [ -f "${coherence_scores}" ]; then
        cp -L "${coherence_scores}" LemonTree/Networks/Module_coherence_scores.txt 2>/dev/null || true
    fi
    
    # Copy enrichment results - skip the large PNG directories to avoid memory issues
    if [ -d "${enrichment_results}" ]; then
        # Only copy CSV files, not PNG images
        find "${enrichment_results}" -maxdepth 3 -name "*.csv" -exec cp --parents {} LemonTree/Enrichment/ \\; 2>/dev/null || true
    elif [ -f "${enrichment_results}" ]; then
        cp -L "${enrichment_results}" LemonTree/Enrichment/ 2>/dev/null || true
    fi
    
    # Copy PKN results
    for f in ${pkn_results}; do
        if [ -d "\$f" ]; then
            cp -rL "\$f"/* LemonTree/PKN_Evaluation/ 2>/dev/null || true
        elif [ -f "\$f" ]; then
            cp -L "\$f" LemonTree/PKN_Evaluation/ 2>/dev/null || true
        fi
    done
    
    # Copy preprocessing results
    for f in ${preprocessing_results}; do
        if [ -d "\$f" ]; then
            cp -rL "\$f"/* LemonTree/Preprocessing/ 2>/dev/null || true
        elif [ -f "\$f" ]; then
            cp -L "\$f" LemonTree/Preprocessing/ 2>/dev/null || true
        fi
    done
    
    # Copy clustering results (Lemon_out directory with tight_clusters and cluster dirs)
    mkdir -p LemonTree/Lemon_out
    for f in ${clustering_results}; do
        if [ -d "\$f" ]; then
            # If it's a directory, copy its contents
            cp -rL "\$f"/* LemonTree/Lemon_out/ 2>/dev/null || true
        elif [ -f "\$f" ]; then
            # If it's a file, copy it directly
            cp -L "\$f" LemonTree/Lemon_out/ 2>/dev/null || true
        fi
    done
    
    # Copy Module_Overview results (Module_Overview.csv, module_expression_analysis.csv, etc.)
    mkdir -p LemonTree/Module_Overview
    for f in ${overview_results}; do
        if [ -d "\$f" ]; then
            # If it's a directory, copy its contents
            cp -rL "\$f"/* LemonTree/Module_Overview/ 2>/dev/null || true
        elif [ -f "\$f" ]; then
            # If it's a file, copy it directly
            cp -L "\$f" LemonTree/Module_Overview/ 2>/dev/null || true
        fi
    done
    
    # Copy parameters log
    if [ -f "${parameters_log}" ]; then
        cp -L "${parameters_log}" ./pipeline_parameters_log.txt 2>/dev/null || true
    fi
    
    # Debug: list what we have
    echo "=== LemonTree directory structure ==="
    find LemonTree -type f | head -50 || true
    echo "======================================"
    
    # Run summary report generator
    if [ -f "${projectDir}/scripts/generate_summary_report.py" ]; then
        SCRIPT_PATH="${projectDir}/scripts/generate_summary_report.py"
    elif [ -f "/app/scripts/generate_summary_report.py" ]; then
        SCRIPT_PATH="/app/scripts/generate_summary_report.py"
    else
        echo "Error: generate_summary_report.py not found"
        exit 1
    fi
    
    # Pass current directory as output_dir with run_id="." so script looks for LemonTree directly
    # The LemonTree directory is already copied to the current working directory by this point
    python3 \$SCRIPT_PATH \\
        --input_dir ${params.input_dir} \\
        --output_dir . \\
        --run_id . \\
        --regulator_types "${regulator_types}" \\
        --parameters_file ./pipeline_parameters_log.txt \\
        --organism ${params.organism}
    
    # Move the report to the expected output location
    if [ -f "${run_id}/Lemonite_Summary_Report.html" ]; then
        mv "${run_id}/Lemonite_Summary_Report.html" ./Lemonite_Summary_Report.html
    fi
    
    echo "Summary report generation completed"
    """

    stub:
    """
    touch Lemonite_Summary_Report.html
    """
}

