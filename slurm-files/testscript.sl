#!/bin/env bash

#SBATCH -J ogledata
#SBATCH -A uoa00357
#SBATCH --time=06:00:00  ## 12 hours max for LPV
#SBATCH --mem-per-cpu=8000
##SBATCH --ntasks=16  ## Use for MPI
##SBATCH --cpus-per-task=1  ## Use for OpenMP
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH -C sb

module load Python/2.7.9-intel-2015a
source ~/mypython/bin/activate

srun echo "running on : $(hostname)"
srun echo "allocation : $SLURM_NODELIST"
srun echo "start time : $(date)"

srun python testscript.py

srun echo "end time : $(date)"
