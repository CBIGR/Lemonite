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
    path "LemonTree/Preprocessing/DESeq_groups.txt", emit: metadata
    path "LemonTree/Preprocessing/name_mapping.tsv", emit: name_mapping
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

    # Determine which preprocessing script to use
    # "proteomics" mode: uses Preprocessing_TFA_Proteomics.R (pre-scaled, Pareto/z-score, TFA optional)
    # "rna" mode (default): uses Preprocessing_TFA_RNA.R (DESeq2 normalisation, TFA enabled)
    if [ "${params.preprocessing_type}" = "proteomics" ] && [ -f "${projectDir}/scripts/Preprocessing_TFA_Proteomics.R" ]; then
        SCRIPT_PATH="${projectDir}/scripts/Preprocessing_TFA_Proteomics.R"
        echo "Using Preprocessing_TFA_Proteomics.R (pre-scaled proteomics mode)"
        
        # Build optional file arguments
        OPTIONAL_ARGS=""
        
        # Check for ID mapping file
        if [ -f "data/Proteins_ID_mapping.tsv" ]; then
            OPTIONAL_ARGS="\$OPTIONAL_ARGS --id_mapping data/Proteins_ID_mapping.tsv"
        fi
        
        # Parse regulator_types for hPTM and metabolomics files
        # Format: "hPTMs:Histone_proteomics.csv,Metabolites:metabolomics_combined.csv"
        REGULATOR_TYPES="${params.regulator_types ?: ''}"
        if [[ "\$REGULATOR_TYPES" == *"hPTM"* ]]; then
            HPTM_FILE=\$(echo "\$REGULATOR_TYPES" | grep -oP '(?<=hPTMs:)[^,]+')
            if [ -f "data/\$HPTM_FILE" ]; then
                OPTIONAL_ARGS="\$OPTIONAL_ARGS --hptm_file data/\$HPTM_FILE"
            fi
        fi
        if [[ "\$REGULATOR_TYPES" == *"Metabolites"* ]]; then
            METABO_FILE=\$(echo "\$REGULATOR_TYPES" | grep -oP '(?<=Metabolites:)[^,:]+')
            if [ -f "data/\$METABO_FILE" ]; then
                OPTIONAL_ARGS="\$OPTIONAL_ARGS --metabolomics_file data/\$METABO_FILE"
            fi
            # Check for metabolomics labels file (explicit param takes priority, then data/ fallback)
            if [ -n "${params.metabolomics_labels_file ?: ''}" ] && [ -f "${params.metabolomics_labels_file}" ]; then
                OPTIONAL_ARGS="\$OPTIONAL_ARGS --metabolomics_labels '${params.metabolomics_labels_file}'"
            elif [ -f "data/metabolomics_name_map.csv" ]; then
                OPTIONAL_ARGS="\$OPTIONAL_ARGS --metabolomics_labels data/metabolomics_name_map.csv"
            fi
        fi
        
        Rscript \$SCRIPT_PATH \\
            --expression "\$EXPRESSION_FILE" \\
            --metadata "\$METADATA_FILE" \\
            --output_dir . \\
            --top_n_genes "${params.top_n_genes}" \\
            --perform_TFA "${params.perform_tfa}" \\
            --organism "${params.organism}" \\
            --sample_id_col "${params.sample_id_col}" \\
            \$OPTIONAL_ARGS
    else
        # Default: RNA-seq preprocessing with DESeq2 and TFA
        if [ -f "${projectDir}/scripts/Preprocessing_TFA_RNA.R" ]; then
            SCRIPT_PATH="${projectDir}/scripts/Preprocessing_TFA_RNA.R"
        elif [ -f "/app/scripts/Preprocessing_TFA_RNA.R" ]; then
            SCRIPT_PATH="/app/scripts/Preprocessing_TFA_RNA.R"
        else
            echo "Error: Neither preprocessing script found (Preprocessing_TFA_RNA.R)"
            exit 1
        fi
        echo "Using Preprocessing_TFA_RNA.R for RNA-seq count data"
        
        # Export gene annotation file path if provided
        export GENE_ANNOTATION_FILE="${params.gene_annotation_file ?: ''}"
        export GENE_ANNOTATION_BACKUP="${projectDir}/PKN/ensembl_mapping_jan2024.txt"

        Rscript \$SCRIPT_PATH \\
            --expression "\$EXPRESSION_FILE" \\
            --metadata "\$METADATA_FILE" \\
            --output_dir . \\
            --regulator_types "${params.regulator_types}" \\
            --top_n_genes "${params.top_n_genes}" \
            --perform_TFA "${params.perform_tfa}" \
            --use_omics_specific_scaling "${params.use_omics_specific_scaling}" \
            --DESeq_contrast1 "${params.deseq_contrast1}" \
            --design_formula "${params.design_formula}" \
            --metadata_columns "${params.metadata_columns}" \
            --expression_col "${params.expression_col}" \
            --sample_id_col "${params.sample_id_col}" \
            --organism "${params.organism}" \
            \$([ ! -z "\$PRIOR_NETWORK" ] && echo "--prior_network \$PRIOR_NETWORK")
    fi
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
    touch LemonTree/Preprocessing/name_mapping.tsv
    """
}
