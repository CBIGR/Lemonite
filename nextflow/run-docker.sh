#!/bin/bash

# Run script for Lemonite Docker pipeline

set -e

# Default values
DATA_DIR=""
OUTPUT_DIR="./results"
IMAGE_NAME="lemontree-pipeline:latest"
INTERACTIVE=false
DOCKER_COMPOSE=false

# Function to show usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Run the Lemonite Analysis Pipeline using Docker"
    echo ""
    echo "OPTIONS:"
    echo "  -d, --data-dir DIR       Input data directory (required)"
    echo "  -o, --output-dir DIR     Output directory (default: ./results)"
    echo "  -i, --interactive        Run in interactive mode"
    echo "  -c, --compose            Use docker-compose"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "EXAMPLES:"
    echo "  # Basic run:"
    echo "  $0 -d /path/to/data"
    echo ""
    echo "  # Custom output directory:"
    echo "  $0 -d /path/to/data -o /path/to/output"
    echo ""
    echo "  # Interactive mode:"
    echo "  $0 -d /path/to/data -i"
    echo ""
    echo "  # Using docker-compose:"
    echo "  $0 -d /path/to/data -c"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -i|--interactive)
            INTERACTIVE=true
            shift
            ;;
        -c|--compose)
            DOCKER_COMPOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$DATA_DIR" ]; then
    echo "Error: Data directory is required"
    usage
    exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory does not exist: $DATA_DIR"
    exit 1
fi

# Convert to absolute paths
DATA_DIR=$(realpath "$DATA_DIR")
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "LemonTree Analysis Pipeline - Docker Runner"
echo "=========================================="
echo "Data directory: $DATA_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Interactive mode: $INTERACTIVE"
echo "Docker compose: $DOCKER_COMPOSE"
echo ""

if [ "$DOCKER_COMPOSE" = true ]; then
    # Use docker-compose
    echo "Starting pipeline with docker-compose..."
    
    # Update docker-compose.yml with actual paths
    export DATA_DIR_HOST="$DATA_DIR"
    export OUTPUT_DIR_HOST="$OUTPUT_DIR"
    
    if [ "$INTERACTIVE" = true ]; then
        docker-compose run --rm lemontree-dev
    else
        docker-compose up --rm lemontree-pipeline
    fi
    
else
    # Use direct docker run
    DOCKER_CMD="docker run --rm"
    
    # Add volume mounts
    DOCKER_CMD="$DOCKER_CMD -v $DATA_DIR:/data:ro"
    DOCKER_CMD="$DOCKER_CMD -v $OUTPUT_DIR:/results"
    
    # Add working directory
    DOCKER_CMD="$DOCKER_CMD -w /app"
    
    if [ "$INTERACTIVE" = true ]; then
        # Interactive mode
        DOCKER_CMD="$DOCKER_CMD -it $IMAGE_NAME bash"
        echo "Starting interactive Docker container..."
        echo "Inside the container, you can run:"
        echo "  nextflow run main.nf --input_dir /data --output_dir /results -profile docker"
    else
        # Pipeline execution mode
        DOCKER_CMD="$DOCKER_CMD $IMAGE_NAME nextflow run main.nf --input_dir /data --output_dir /results -profile docker"
        echo "Starting pipeline execution..."
    fi
    
    echo "Docker command: $DOCKER_CMD"
    echo ""
    
    # Execute the command
    eval $DOCKER_CMD
fi

echo ""
echo "Pipeline execution completed!"
echo "Results can be found in: $OUTPUT_DIR"
