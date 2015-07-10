#!/bin/bash -l
#SBATCH --job-name=OMP-6
#SBATCH --account=director100
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --export=NONE

export OMP_NUM_THREADS=6

aprun -n 4 -S 2 -d 6 ./PoreProcess 

