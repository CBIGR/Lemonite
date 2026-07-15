#!/bin/bash -l
#PBS -N pkn_pipeline
#PBS -l nodes=1:ppn=24
#PBS -l mem=128gb
#PBS -l walltime=72:00:00
#PBS -m abe
# ============================================================================
# HPC submission script for the LemonIte PKN build pipeline (VSC Ghent / Doduo).
# Submit with:  qsub submit_pkn_hpc.sh
#
# Resumable: if the job hits the walltime, simply re-submit — --resume continues
# from the last checkpointed database without re-querying completed sources.
# ============================================================================

LEMONITE_DIR=/user/gent/435/vsc43501/25BVM_Lemonite
PIPELINE_DIR=$LEMONITE_DIR/repo/LemonIte/build_PKN/PKN_build_pipeline_opus4.8
VENV=$LEMONITE_DIR/venv_pkn_doduo

source "$VENV/bin/activate"

# Pipeline configuration (overrides config.py defaults; see config_hpc.py)
export PKN_WORKDIR=$LEMONITE_DIR
export PKN_OUTPUT_DIR_NAME=PKN
export PKN_DB_DIR=$LEMONITE_DIR/databases
export PKN_GEM_DIR=$LEMONITE_DIR/models/Human1-GEM/model
export PKN_CONFIG=config_hpc

cd "$PIPELINE_DIR"

echo "=============================="
echo "PKN Pipeline - HPC Run"
echo "Started: $(date)"
echo "Node: $HOSTNAME"
echo "Output: $PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME"
echo "=============================="

python main.py --all --resume --workers 24

echo "=============================="
echo "Finished: $(date)"
echo "=============================="
