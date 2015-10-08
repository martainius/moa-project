#!/bin/env bash

#SBATCH -J classifier
#SBATCH -A uoa00357
#SBATCH --time=01:00:00
#SBATCH --mem=36GB
#SBATCH --nodes=1
##SBATCH --cpus-per-task=1
##SBATCH --ntasks-per-node=1
#SBATCH -C sb

##module load Python/2.7.9-intel-2015a
##source ~/mypython/bin/activate

srun echo "running on : $(hostname)"
srun echo "allocation : $SLURM_NODELIST"
srun echo "start time : $(date)"

srun java -Xmx32G -cp /home/mdon849/weka-3-7-11/weka.jar weka.classifiers.trees.RandomForest -I 1000 -t /projects/uoa00357/moa/training/training.arff -d /projects/uoa00357/moa/training/ftrs_2.model > /projects/uoa00357/moa/training/ftrs_2.dat

srun echo "end time : $(date)"
