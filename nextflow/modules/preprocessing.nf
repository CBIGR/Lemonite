process PREPROCESSING_TFA {
    tag "Preprocessing and TFA analysis"

    input:
    path input_dir
    path data_dir
    val run_id

    output:
    path "LemonTree/Preprocessing/LemonPreprocessed_expression.txt", emit: preprocessed_data
    path "LemonTree/Preprocessing/LemonPreprocessed_complete.txt", emit: complete_data
    path "LemonTree/Preprocessing/LemonPreprocessed_*.txt", emit: omics_preprocessed
    path "LemonTree/Preprocessing/*.txt", emit: regulator_files
    path "LemonTree/Preprocessing/DESeq_groups.txt", emit: metadata
    path "LemonTree/Preprocessing/PCA_*.pdf", emit: pca_plots, optional: true
    path "LemonTree/Preprocessing/*", emit: lemontree_inputs
    path "LemonTree/Preprocessing/*", emit: preprocessing_results, optional: true
    path "TFA/*", emit: tfa_results, optional: true

    script:
    """
    # Resolve expression file: explicit param or auto-detect with uniqueness check
    if [ -n "${params.expression_file ?: ''}" ]; then
        EXPRESSION_FILE="${params.expression_file}"
    else
        MATCHES=\$(find data/ -name "*host_tx_counts.tsv" -o -name "*expression*.tsv" -o -name "*counts*.tsv" -o -name "Counts.tsv")
        MATCH_COUNT=\$(echo "\$MATCHES" | grep -c . || true)
        if [ "\$MATCH_COUNT" -eq 0 ]; then
            echo "Error: Expression file not found. Provide --expression_file or place a *counts*.tsv / Counts.tsv in data/"
            find data/ -type f -name "*.tsv" -o -name "*.txt"
            exit 1
        elif [ "\$MATCH_COUNT" -gt 1 ]; then
            echo "Error: Multiple expression files matched (\$MATCH_COUNT). Please specify --expression_file explicitly:"
            echo "\$MATCHES"
            exit 1
        fi
        EXPRESSION_FILE=\$(echo "\$MATCHES" | head -1)
    fi

    # Resolve metadata file: explicit param or auto-detect with uniqueness check
    if [ -n "${params.metadata_file ?: ''}" ]; then
        METADATA_FILE="${params.metadata_file}"
    else
        MATCHES=\$(find data/ -name "*metadata*.txt" -o -name "*Metadata*.txt")
        MATCH_COUNT=\$(echo "\$MATCHES" | grep -c . || true)
        if [ "\$MATCH_COUNT" -eq 0 ]; then
            echo "Error: Metadata file not found. Provide --metadata_file or place a *metadata*.txt in data/"
            exit 1
        elif [ "\$MATCH_COUNT" -gt 1 ]; then
            echo "Error: Multiple metadata files matched (\$MATCH_COUNT). Please specify --metadata_file explicitly:"
            echo "\$MATCHES"
            exit 1
        fi
        METADATA_FILE=\$(echo "\$MATCHES" | head -1)
    fi

    # Resolve prior network (optional)
    PRIOR_NETWORK=\$(find data/ -name "*CollecTRI*.txt" -o -name "*network*.txt" | head -1)
    if [ -z "\$PRIOR_NETWORK" ] && [ -f "${projectDir}/PKN/CollecTRI_network.txt" ]; then
        PRIOR_NETWORK="${projectDir}/PKN/CollecTRI_network.txt"
        echo "Using bundled network from projectDir: \$PRIOR_NETWORK"
    fi

    echo "Found core files:"
    echo "Expression: \$EXPRESSION_FILE"
    echo "Metadata: \$METADATA_FILE"
    echo "Prior network: \$PRIOR_NETWORK"

    # Copy name_map.csv if it exists
    if [ -f "data/name_map.csv" ]; then
        mkdir -p "LemonTree/Preprocessing/"
        cp "data/name_map.csv" "LemonTree/Preprocessing/"
        echo "Copied name_map.csv to Preprocessing directory"
    fi

    # Run R preprocessing script with regulator_types parameter
    # Check host script first (allows updates without rebuilding container)
    if [ -f "${projectDir}/scripts/preprocessing_tfa_complete.R" ]; then
        SCRIPT_PATH="${projectDir}/scripts/preprocessing_tfa_complete.R"
    elif [ -f "/app/scripts/preprocessing_tfa_complete.R" ]; then
        SCRIPT_PATH="/app/scripts/preprocessing_tfa_complete.R"
    else
        echo "Error: preprocessing_tfa_complete.R not found"
        exit 1
    fi
    Rscript \$SCRIPT_PATH \\
        --expression "\$EXPRESSION_FILE" \\
        --metadata "\$METADATA_FILE" \\
        --output_dir . \\
        --regulator_types "${params.regulator_types}" \\
        --top_n_genes ${params.top_n_genes} \\
        --perform_TFA ${params.perform_tfa} \\
        --use_omics_specific_scaling ${params.use_omics_specific_scaling} \\
        --DESeq_contrast1 ${params.deseq_contrast1} \\
        --design_formula "${params.design_formula}" \\
        --metadata_columns "${params.metadata_columns}" \\
        --expression_col ${params.expression_col} \\
        --sample_id_col "${params.sample_id_col}" \\
        --organism ${params.organism} \\
        \$([ ! -z "\$PRIOR_NETWORK" ] && echo "--prior_network \$PRIOR_NETWORK")
    """

    stub:
    """
    mkdir -p LemonTree/Preprocessing TFA
    touch LemonTree/Preprocessing/LemonPreprocessed_expression.txt
    touch LemonTree/Preprocessing/LemonPreprocessed_complete.txt
    touch LemonTree/Preprocessing/LemonPreprocessed_metabolomics.txt
    touch LemonTree/Preprocessing/DESeq_groups.txt
    touch LemonTree/Preprocessing/lovering_TFs.txt
    touch LemonTree/Preprocessing/metabolites.txt
    """
}
