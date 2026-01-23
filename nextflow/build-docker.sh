#!/bin/bash

# Build script for LemonTree Docker image

set -e

echo "Building Lemonite Analysis Pipeline Docker image..."
echo "=================================================="

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Build arguments
IMAGE_NAME="lemontree-pipeline"
IMAGE_TAG="latest"
FULL_IMAGE_NAME="${IMAGE_NAME}:${IMAGE_TAG}"

# Build the Docker image
echo "Building Docker image: $FULL_IMAGE_NAME"
docker build -t "$FULL_IMAGE_NAME" .

# Check if build was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Docker image built successfully!"
    echo ""
    echo "Image details:"
    docker images "$IMAGE_NAME" --format "table {{.Repository}}\\t{{.Tag}}\\t{{.ID}}\\t{{.CreatedAt}}\\t{{.Size}}"
    echo ""
    echo "Usage examples:"
    echo "  # Run pipeline:"
    echo "  docker run -v /path/to/data:/data -v ./results:/results $FULL_IMAGE_NAME nextflow run main.nf --input_dir /data --output_dir /results -profile docker"
    echo ""
    echo "  # Interactive shell:"
    echo "  docker run -it -v /path/to/data:/data $FULL_IMAGE_NAME bash"
    echo ""
    echo "  # Using docker-compose:"
    echo "  docker-compose up lemontree-pipeline"
    echo ""
else
    echo "❌ Docker build failed!"
    exit 1
fi
