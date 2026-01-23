#!/bin/bash

# Build script for Lemonite Singularity container

set -e

echo "Building Lemonite Analysis Pipeline Singularity image..."
echo "=================================================="

# Check if singularity is installed
if ! command -v singularity &> /dev/null; then
    echo "Error: Singularity is not installed or not in PATH"
    echo "Please install Singularity first: https://sylabs.io/guides/3.0/user-guide/installation.html"
    exit 1
fi

# Create user cache and temp directories to avoid filling up system root
export SINGULARITY_CACHEDIR="${PWD}/.singularity/cache"
export SINGULARITY_TMPDIR="${PWD}/.singularity/tmp"
export APPTAINER_CACHEDIR="${PWD}/.singularity/cache"
export APPTAINER_TMPDIR="${PWD}/.singularity/tmp"
export TMPDIR="${PWD}/.singularity/tmp"
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

echo "Using cache directory: $SINGULARITY_CACHEDIR"
echo "Using temp directory: $SINGULARITY_TMPDIR"

# Check available space in user directory
USER_SPACE=$(df -h "$PWD" | awk 'NR==2 {print $4}')
echo "Available space in user directory: $USER_SPACE"

# Build the image
echo "Building Singularity image: lemontree-pipeline.sif"
echo "Note: Building without sudo using --fakeroot (requires user namespace support)"
singularity build --fakeroot --bind "${PWD}/.singularity/tmp:/tmp" lemontree-pipeline.sif Singularity.def

# Cleanup function
cleanup_cache() {
    echo "Cleaning up build cache..."
    if [ -d "$SINGULARITY_CACHEDIR" ]; then
        du -sh "$SINGULARITY_CACHEDIR" 2>/dev/null || echo "Cache directory size: unknown"
        rm -rf "$SINGULARITY_CACHEDIR"
        echo "Cache directory removed"
    fi
}

# Verify the build
if [ -f "lemontree-pipeline.sif" ]; then
    echo ""
    echo "✅ Singularity image built successfully!"
    echo "Image size: $(du -h lemontree-pipeline.sif | cut -f1)"
    echo ""
    cleanup_cache
    echo ""
    echo "Test the image with:"
    echo "singularity exec lemontree-pipeline.sif R --version"
    echo "singularity exec lemontree-pipeline.sif python3 --version"
    echo ""
    echo "Run the pipeline with:"
    echo "nextflow run main.nf -profile singularity --input_dir /path/to/data"
else
    echo "❌ Failed to build Singularity image"
    exit 1
fi
