#!/bin/bash -l
#SBATCH --job-name=OMP-2
#SBATCH --account=director100
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --export=NONE

export OMP_NUM_THREADS=2

aprun -n 12 -S 6 -d 2 ./PoreProcess 

