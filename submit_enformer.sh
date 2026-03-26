#!/bin/bash

#SBATCH -p gpu
#SBATCH --gres=gpu:A6000:1          
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name enformer
#SBATCH --output /home/azstephe/liverRegression/regression_liver/data/slurm_out/enformer-%J.o
#SBATCH --error=/home/azstephe/liverRegression/regression_liver/data/slurm_out/enformer-%j.e
#SBATCH -t 0-08:00:00

source /home/azstephe/miniconda3/etc/profile.d/conda.sh
source activate enformer_env

mkdir -p /home/azstephe/liverRegression/regression_liver/data/enformer_outputs

python /home/azstephe/liverRegression/regression_liver/scripts/eval_enformer.py
