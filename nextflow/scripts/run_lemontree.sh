#!/bin/bash

# LemonTree clustering script for Nextflow
# Adapted from lemon_cluster.pbs

# Parameters
CLUSTER_ID=$1
INPUT_FILE=$2
OUTPUT_DIR=$3
LEMONTREE_JAR=$4

echo "=== LemonTree Clustering Debug Info ==="
echo "Cluster ID: ${CLUSTER_ID}"
echo "Input file: ${INPUT_FILE}"
echo "Output dir: ${OUTPUT_DIR}"
echo "LemonTree JAR: ${LEMONTREE_JAR}"
echo "========================================"

# Create output directory
mkdir -p "${OUTPUT_DIR}/Lemon_results"

# Run LemonTree clustering
if command -v java &> /dev/null; then
    echo "Java found: $(java -version 2>&1 | head -n 1)"
    
    # Build classpath including all dependencies
    LEMONTREE_DIR=$(dirname "${LEMONTREE_JAR}")
    echo "LemonTree directory: ${LEMONTREE_DIR}"
    
    CLASSPATH="${LEMONTREE_JAR}"
    if [ -d "${LEMONTREE_DIR}/lib" ]; then
        echo "Found lib directory at ${LEMONTREE_DIR}/lib"
        for jar in "${LEMONTREE_DIR}"/lib/*.jar; do
            CLASSPATH="${CLASSPATH}:${jar}"
        done
        echo "Classpath includes $(ls -1 "${LEMONTREE_DIR}"/lib/*.jar 2>/dev/null | wc -l) dependency JARs"
    else
        echo "Warning: lib directory not found at ${LEMONTREE_DIR}/lib"
    fi
    
    echo "Final classpath: ${CLASSPATH}"
    echo "Running LemonTree..."
    
    # LemonTree clustering - save to Lemon_results subdirectory
    java -cp "${CLASSPATH}" lemontree.modulenetwork.RunCli -task ganesh -data_file "${INPUT_FILE}" -output_file "${OUTPUT_DIR}/Lemon_results/cluster_${CLUSTER_ID}"
    
    JAVA_EXIT_CODE=$?
    if [ "${JAVA_EXIT_CODE}" -ne 0 ]; then
        echo "Error: LemonTree clustering failed with exit code ${JAVA_EXIT_CODE}"
        exit "${JAVA_EXIT_CODE}"
    fi

else
    echo "Error: Java not found. Please ensure Java is installed and available."
    exit 1
fi

echo "LemonTree clustering completed for cluster ${CLUSTER_ID}"
