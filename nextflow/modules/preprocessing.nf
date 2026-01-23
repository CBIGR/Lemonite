process PREPROCESSING_TFA {
    tag "Preprocessing and TFA analysis"
    publishDir "${(params.output_dir ?: (params.input_dir ? params.input_dir + '/results' : './results'))}/${(params.computed_run_id ?: (params.run_id ?: 'run_auto'))}", mode: 'copy'

    input:
    path input_dir
    path data_dir
    val run_id

    output:
    path "LemonTree/Preprocessing/LemonPreprocessed_expression.txt", emit: preprocessed_data
    path "LemonTree/Preprocessing/LemonPreprocessed_complete.txt", emit: complete_data
    path "LemonTree/Preprocessing/*.txt", emit: regulator_files
    path "LemonTree/Preprocessing/DESeq_groups.txt", emit: metadata
    path "LemonTree/Preprocessing/PCA_*.pdf", emit: pca_plots, optional: true
    path "LemonTree/Preprocessing/*", emit: lemontree_inputs
    path "LemonTree/Preprocessing/*", emit: preprocessing_results, optional: true
    path "TFA/*", emit: tfa_results, optional: true

    script:
    """
    # Find the required input files in the data directory
    EXPRESSION_FILE=\$(find data/ -name "*host_tx_counts.tsv" -o -name "*expression*.tsv" -o -name "*counts*.tsv" -o -name "Counts.tsv" | head -1)
    METADATA_FILE=\$(find data/ -name "*metadata*.txt" -o -name "*Metadata*.txt" | head -1)
    PRIOR_NETWORK=\$(find data/ -name "*CollecTRI*.txt" -o -name "*network*.txt" | head -1)

    # Check if required files exist
    if [ -z "\$EXPRESSION_FILE" ]; then
        echo "Error: Expression file not found"
        echo "Available files:"
        find data/ -type f -name "*.tsv" -o -name "*.txt"
        exit 1
    fi
    if [ -z "\$METADATA_FILE" ]; then
        echo "Error: Metadata file not found in ${data_dir}"
        exit 1
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
    mkdir -p LemonTree Preprocessing TFA
    touch LemonTree/LemonPreprocessed_expression.txt
    touch LemonTree/LemonPreprocessed_complete.txt
    touch LemonTree/LemonPreprocessed_metabolomics.txt
    touch LemonTree/DESeq_groups.txt
    touch LemonTree/lovering_TFs.txt
    touch LemonTree/metabolites.txt
    """
}
