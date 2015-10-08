#!/bin/env bash

#SBATCH -J moadata
#SBATCH -A uoa00357
#SBATCH --time=09:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH -C sb
#SBATCH --array=1-220

module load Python/2.7.9-intel-2015a
source ~/mypython/bin/activate

srun echo "running on : $(hostname)"
srun echo "allocation : $SLURM_NODELIST"
srun echo "start time : $(date)"

PARAMETERS=$(sed -n ${SLURM_ARRAY_TASK_ID}p moa_params.dat)

srun python testData.py $PARAMETERS

srun echo "end time : $(date)"
