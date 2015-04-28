#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --account=director100
#SBATCH --job-name=test
#SBATCH --time=01:00:00
#SBATCH --partition=debugq
#SBATCH --export=NONE

export OMP_NUM_THREADS=12
aprun -n 1 -d 12 ../src/PoreProcess+pat
#aprun -n 1 ../src/PoreProcess

