#!/bin/bash -l
#SBATCH --job-name=OMP-4
#SBATCH --account=director100
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --export=NONE

export OMP_NUM_THREADS=4

aprun -n 6 -S 3 -d 4 ./PoreProcess 

