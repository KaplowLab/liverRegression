#!/bin/bash

#SBATCH -p 
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --job-name enformer
#SBATCH --output 
#SBATCH --error
#SBATCH -t 0-08:00:00

source /home/azstephe/miniconda3/etc/profile.d/conda.sh
source activate enformer_env

mkdir -p /home/azstephe/liverRegression/regression_liver/data/enformer_outputs

python /home/azstephe/liverRegression/regression_liver/scripts/eval_enformer.py
