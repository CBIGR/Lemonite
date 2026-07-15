#!/bin/bash -l
#PBS -N pkn_pipeline_test
#PBS -l nodes=1:ppn=4
#PBS -l mem=64gb
#PBS -l walltime=02:00:00
#PBS -m abe

# ============================================================
# TEST HPC submission script for the PKN Build Pipeline
# Runs Step 1 only with a small metabolite subset and
# fast local-file-based databases (no heavy API calls).
# Submit with: qsub submit_pkn_hpc_test.sh
# ============================================================

LEMONITE_DIR=/user/gent/435/vsc43501/25BVM_Lemonite
PIPELINE_DIR=$LEMONITE_DIR/repo/LemonIte/build_PKN/PKN_build_pipeline
VENV=$LEMONITE_DIR/venv_pkn_doduo

# Activate virtual environment
source $VENV/bin/activate

# Set pipeline configuration via environment variables (overrides config.py defaults)
export PKN_WORKDIR=$LEMONITE_DIR
export PKN_OUTPUT_DIR_NAME=PKN_test
export PKN_DB_DIR=$LEMONITE_DIR/databases
export PKN_GEM_DIR=$LEMONITE_DIR/models/Human1-GEM/model
export PKN_CONFIG=config_hpc

cd $PIPELINE_DIR

echo "=============================="
echo "PKN Pipeline - TEST Run"
echo "Started: $(date)"
echo "Node: $HOSTNAME"
echo "PKN_WORKDIR:  $PKN_WORKDIR"
echo "PKN_DB_DIR:   $PKN_DB_DIR"
echo "PKN_GEM_DIR:  $PKN_GEM_DIR"
echo "Output dir:   $PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME"
echo "=============================="

# Test settings:
#   --step 1            : run only the metabolite-gene interaction step
#   --max-metabolites 50: limit to first 50 metabolites for speed
#   --databases         : only local-file databases (no slow API retrievers)
#   --workers 4         : match available CPUs
python main.py \
    --step 1 \
    --max-metabolites 50 \
    --databases biogrid,stitch,lincs,metalinks,gem_dist1 \
    --workers 4 \
    --output-dir PKN_test

echo "=============================="
echo "Finished: $(date)"
echo "Exit code: $?"
echo "=============================="
