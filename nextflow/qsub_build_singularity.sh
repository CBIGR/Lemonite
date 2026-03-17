#!/bin/bash -l
#PBS -l nodes=1:ppn=1
#PBS -l walltime=04:00:00
#PBS -l mem=16gb
#PBS -m abe

cd $PBS_O_WORKDIR

./build-singularity.sh
