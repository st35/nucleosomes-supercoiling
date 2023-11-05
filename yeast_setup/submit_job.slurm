#!/bin/bash

#SBATCH --job-name=gapr
#SBATCH --ntasks=1
#SBATCH --time=07:59:59

python Get_Avg_Gene_Expression.py

echo $SLURM_JOB_NODELIST
