#!/bin/bash
# Quick start script for the PKN pipeline

echo "========================================"
echo "PKN Pipeline - Quick Start"
echo "========================================"

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install dependencies
echo "Installing dependencies..."
pip install -r requirements.txt

echo ""
echo "Setup complete!"
echo ""
echo "Available commands:"
echo "  python main.py --all                    # Run full pipeline"
echo "  python main.py --step 1                 # Run Step 1 only"
echo "  python main.py --step 1 --databases biogrid,stitch  # Specific DBs"
echo "  python main.py --step 1 --resume        # Resume from checkpoint"
echo ""
echo "Test individual retriever:"
echo "  python step1_metabolites/biogrid.py"
echo ""
