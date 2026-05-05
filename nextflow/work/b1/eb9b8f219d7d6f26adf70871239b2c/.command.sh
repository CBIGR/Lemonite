#!/bin/bash -euo pipefail
# Validate input file
if [ ! -s "LemonPreprocessed_expression.txt" ]; then
    echo "Error: Input file 'LemonPreprocessed_expression.txt' is missing or empty"
    exit 1
fi

# Run LemonTree clustering using our script
# Check host script first (allows updates without rebuilding container)
if [ -f "/home/borisvdm/repo/LemonIte/nextflow/scripts/run_lemontree.sh" ]; then
    SCRIPT_PATH="/home/borisvdm/repo/LemonIte/nextflow/scripts/run_lemontree.sh"
elif [ -f "/app/scripts/run_lemontree.sh" ]; then
    SCRIPT_PATH="/app/scripts/run_lemontree.sh"
else
    echo "Error: run_lemontree.sh not found"
    exit 1
fi
bash $SCRIPT_PATH         3         LemonPreprocessed_expression.txt         .         "/home/borisvdm/repo/LemonIte/nextflow/LemonTree/lemontree_v3.1.1.jar"         42
