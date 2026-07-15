#!/bin/bash -l
#PBS -N pkn_test_500
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -l walltime=04:00:00
#PBS -m abe
# ============================================================================
# HPC SMOKE TEST for the LemonIte PKN build pipeline (VSC Ghent / Doduo).
#
# Runs step 1 on the first 500 metabolites into a throwaway output dir
# (PKN_test_500) — the validation build described in HANDOFF.md. Use this to
# confirm the environment, database paths and API access work before
# committing to the full multi-day run (submit_pkn_hpc.sh).
#
# Submit with:  qsub submit_pkn_test.sh
# ============================================================================

LEMONITE_DIR=/user/gent/435/vsc43501/25BVM_Lemonite
PIPELINE_DIR=$LEMONITE_DIR/repo/LemonIte/build_PKN/PKN_build_pipeline_opus4.8
VENV=$LEMONITE_DIR/venv_pkn_doduo

source "$VENV/bin/activate"

# Pipeline configuration (overrides config.py defaults; see config_hpc.py)
export PKN_WORKDIR=$LEMONITE_DIR
export PKN_OUTPUT_DIR_NAME=PKN_test_500
export PKN_DB_DIR=$LEMONITE_DIR/databases
export PKN_GEM_DIR=$LEMONITE_DIR/models/Human1-GEM/model
export PKN_CONFIG=config_hpc

cd "$PIPELINE_DIR"

echo "=============================="
echo "PKN Pipeline - HPC SMOKE TEST (step 1, 500 metabolites)"
echo "Started: $(date)"
echo "Node: $HOSTNAME"
echo "Output: $PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME"
echo "=============================="

python main.py --step 1 --max-metabolites 500 --output-dir PKN_test_500 --workers 8

echo "=============================="
echo "Finished: $(date)"
echo "=============================="
