#!/bin/bash -l
#SBATCH --job-name=OMP-3
#SBATCH --account=director100
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --export=NONE

export OMP_NUM_THREADS=3

aprun -n 8 -S 4 -d 3 ./PoreProcess 

