#!/bin/bash -l
# ============================================================================
# Idempotent monthly PKN rebuild wrapper (for cron or manual invocation).
#
# Safe to re-trigger: --resume skips already-completed databases, so a second
# run shortly after a successful one is a fast no-op, and a run after an
# interrupted/throttled attempt finishes where it left off.
#
# Configure the paths below (or rely on the environment already exporting the
# PKN_* variables), then schedule via schedule/crontab.example.
# ============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PIPELINE_DIR"

# --- Configuration (override via environment if already set) ----------------
export PKN_WORKDIR="${PKN_WORKDIR:-$HOME/Lemonite}"
export PKN_OUTPUT_DIR_NAME="${PKN_OUTPUT_DIR_NAME:-PKN}"
export PKN_DB_DIR="${PKN_DB_DIR:-$PKN_WORKDIR/databases}"
export PKN_GEM_DIR="${PKN_GEM_DIR:-$PKN_WORKDIR/models/Human1-GEM/model}"

# Activate the virtual environment if present.
if [ -f "$PIPELINE_DIR/venv/bin/activate" ]; then
    # shellcheck disable=SC1091
    source "$PIPELINE_DIR/venv/bin/activate"
fi

STAMP="$(date +%Y%m%d)"
echo "=============================================="
echo "PKN monthly rebuild — $(date)"
echo "  PKN_WORKDIR=$PKN_WORKDIR"
echo "  PKN_DB_DIR=$PKN_DB_DIR"
echo "  output: $PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME"
echo "=============================================="

python main.py --all --resume

echo "Finished monthly rebuild ($STAMP): $(date)"
