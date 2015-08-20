#!/bin/bash -l
#SBATCH --job-name=Job1
#SBATCH --account=director100
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --partition=debugq
#SBATCH --export=NONE

export OMP_NUM_THREADS=1

aprun -n 4 ./PoreProcess 

