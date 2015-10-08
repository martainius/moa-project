#!/bin/env bash

#SBATCH -J classifier
#SBATCH -A uoa00357
#SBATCH --time=00:30:00
#SBATCH --mem=48GB
##SBATCH --mem-per-cpu=1000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
##SBATCH --ntasks-per-node=1
#SBATCH -C sb

##module load Python/2.7.9-intel-2015a
##source ~/mypython/bin/activate

srun echo "running on : $(hostname)"
srun echo "allocation : $SLURM_NODELIST"
srun echo "start time : $(date)"

srun java -Xmx44G -cp ~/weka-3-7-11/weka.jar weka.core.SystemInfo

srun echo "end time : $(date)"
