process CLUSTERING {
    tag "CLUSTERING $cluster_id"
    
    input:
    tuple path(preprocessed_data), val(cluster_id)
    val run_id

    output:
    path "Lemon_results/cluster_${cluster_id}", emit: clusters

    script:
    """
    # Run LemonTree clustering using our script
    # Check host script first (allows updates without rebuilding container)
    if [ -f "${projectDir}/scripts/run_lemontree.sh" ]; then
        SCRIPT_PATH="${projectDir}/scripts/run_lemontree.sh"
    elif [ -f "/app/scripts/run_lemontree.sh" ]; then
        SCRIPT_PATH="/app/scripts/run_lemontree.sh"
    else
        echo "Error: run_lemontree.sh not found"
        exit 1
    fi
    bash \$SCRIPT_PATH \
        ${cluster_id} \
        ${preprocessed_data} \
        . \
        ${params.lemontree_jar}
    """

    stub:
    """
    mkdir -p Lemon_results/cluster_${cluster_id}
    touch Lemon_results/cluster_${cluster_id}/result.txt
    """
}
