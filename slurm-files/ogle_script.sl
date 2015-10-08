#!/bin/env bash

#SBATCH -J ogledata
#SBATCH -A uoa00357
#SBATCH --time=08:50:00  ## 12 hours max for LPV
#SBATCH --mem-per-cpu=3000
##SBATCH --ntasks=16  ## Use for MPI
##SBATCH --cpus-per-task=1  ## Use for OpenMP
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH -C sb
#SBATCH --array=1-16  ## The size of this job array should be equivalent to 
                      ## $CORES in split.sh (i.e. the number of CPUs each csv 
                      ## file is split over).

module load Python/2.7.9-intel-2015a
source ~/mypython/bin/activate

srun echo "running on : $(hostname)"
srun echo "allocation : $SLURM_NODELIST"
srun echo "start time : $(date)"

srun python synthData_long.py $1-$SLURM_ARRAY_TASK_ID

srun echo "end time : $(date)"
