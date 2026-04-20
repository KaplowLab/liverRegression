#!/bin/bash

#SBATCH -p 
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --job-name=alpha
#SBATCH --output 
#SBATCH --error
#SBATCH -t 0-12:00:00

conda init bash
source ~/.bashrc
conda activate alphagenome_env

python /home/azstephe/liverRegression/regression_liver/scripts/run_alphagenome.py
