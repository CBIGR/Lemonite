#!/bin/bash -l
#PBS -N pkn_pipeline
#PBS -l nodes=1:ppn=24
#PBS -l mem=128gb
#PBS -l walltime=72:00:00
#PBS -m abe

# ============================================================
# HPC submission script for the PKN Build Pipeline
# Submit with: qsub submit_pkn_hpc.sh
# ============================================================

LEMONITE_DIR=/user/gent/435/vsc43501/25BVM_Lemonite
PIPELINE_DIR=$LEMONITE_DIR/repo/LemonIte/build_PKN/PKN_build_pipeline
VENV=$LEMONITE_DIR/venv_pkn_doduo

# Activate virtual environment
source $VENV/bin/activate

# Set pipeline configuration via environment variables (overrides config.py defaults)
export PKN_WORKDIR=$LEMONITE_DIR
export PKN_OUTPUT_DIR_NAME=PKN
export PKN_DB_DIR=$LEMONITE_DIR/databases
export PKN_GEM_DIR=$LEMONITE_DIR/models/Human1-GEM/model
export PKN_CONFIG=config_hpc

cd $PIPELINE_DIR

echo "=============================="
echo "PKN Pipeline - HPC Run"
echo "Started: $(date)"
echo "Node: $HOSTNAME"
echo "PKN_WORKDIR:  $PKN_WORKDIR"
echo "PKN_DB_DIR:   $PKN_DB_DIR"
echo "PKN_GEM_DIR:  $PKN_GEM_DIR"
echo "Output dir:   $PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME"
echo "=============================="

python main.py --all --workers 24

echo "=============================="
echo "Finished: $(date)"
echo "=============================="
