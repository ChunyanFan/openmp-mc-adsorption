#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --account=director100
#SBATCH --job-name=test
#SBATCH --time=00:14:00
#SBATCH --partition=workq
#SBATCH --export=NONE

#module swap PrgEnv-intel PrgEnv-cray
#module load perftools
#export LD_LIBRARY_PATH=/group/partner1032/ddeeptimahanti/src:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=10
#export OMP_DYNAMIC=TRUE
#export OMP_WAIT_POLICY=PASSIVE
#export OMP_NESTED=true
#export OMP_SCHEDULE=RUNTIME
#aprun -n 1 ../src/PoreProcess+pat
aprun -n 1 -d 10 -cc 0-9 ../src/PoreProcess
#aprun -n 1 -N 1 -d 10 -cc 0-9 PoreProcess+apa

#aprun -n 1 ../src/PoreProcess
