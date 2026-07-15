#!/bin/bash -l
# Create a virtual environment and install dependencies for the PKN pipeline.
set -euo pipefail

cd "$(dirname "$0")"

VENV="${1:-venv}"

python3 -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"
pip install --upgrade pip
pip install -r requirements.txt

echo
echo "Done. Activate with:  source $(pwd)/$VENV/bin/activate"
echo "Then run:             python main.py --all"
