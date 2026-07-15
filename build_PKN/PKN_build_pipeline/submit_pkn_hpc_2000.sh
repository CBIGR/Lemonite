#!/bin/bash -l
#PBS -N pkn_test_2000
#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb
#PBS -l walltime=12:00:00
#PBS -m abe

# ============================================================
# Mid-scale TEST: first 2000 metabolites, all 3 steps.
# Uses local-file databases only (no slow API calls).
# Submit with: qsub submit_pkn_hpc_2000.sh
# ============================================================

LEMONITE_DIR=/user/gent/435/vsc43501/25BVM_Lemonite
PIPELINE_DIR=$LEMONITE_DIR/repo/LemonIte/build_PKN/PKN_build_pipeline
VENV=$LEMONITE_DIR/venv_pkn_doduo

# Activate virtual environment
source $VENV/bin/activate

# Pipeline configuration
export PKN_WORKDIR=$LEMONITE_DIR
export PKN_OUTPUT_DIR_NAME=PKN_test_2000
export PKN_DB_DIR=$LEMONITE_DIR/databases
export PKN_GEM_DIR=$LEMONITE_DIR/models/Human1-GEM/model
export PKN_CONFIG=config_hpc

cd $PIPELINE_DIR

echo "=============================="
echo "PKN Pipeline - 2000-metabolite test (all steps)"
echo "Started: $(date)"
echo "Node: $HOSTNAME"
echo "PKN_WORKDIR:  $PKN_WORKDIR"
echo "PKN_DB_DIR:   $PKN_DB_DIR"
echo "PKN_GEM_DIR:  $PKN_GEM_DIR"
echo "Output dir:   $PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME"
echo "=============================="

python main.py \
    --all \
    --max-metabolites 2000 \
    --databases biogrid,stitch,lincs,metalinks,gem_dist1 \
    --workers 8 \
    --output-dir PKN_test_2000

echo "=============================="
echo "Finished: $(date)"
echo "Exit code: $?"
