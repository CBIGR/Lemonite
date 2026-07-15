#!/bin/bash -l
#PBS -l nodes=1:ppn=8
#PBS -l walltime=02:00:00
#PBS -l mem=32gb
#PBS -m abe
#PBS -M boris.vandemoortele@ugent.be

cd $PBS_O_WORKDIR

./build-singularity_parallel.sh
