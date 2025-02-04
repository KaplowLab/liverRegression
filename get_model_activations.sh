#!/bin/bash

set -e

cd /home/azstephe/repos/cnn_pipeline

runname=$1
modelname=$2

mkdir -p "/home/azstephe/regression_liver/data/model_outputs/${modelname}"

srun -p pfen3 -n 1 --gres gpu:1 --pty bash -c "

	source ~/.bashrc && \
	conda activate keras2-tf27 && \

	python scripts/get_activations.py \
    	  -model "/home/azstephe/repos/cnn_pipeline/wandb/${runname}/files/model-best.h5" \
    	  -in_files "/home/azstephe/regression_liver/data/splits/positives/mouse_liver_TRAIN.narrowPeak" \
    	  "/home/azstephe/regression_liver/data/splits/negatives/nonMouse_liver_andRat_andCow_andPig_andMacaque_TRAIN_500bp.bed" \
    	  -in_genomes "/projects/pfenninggroup/machineLearningForComputationalBiology/halLiftover_chains/data/raw_data/fasta/Mus_musculus.fa" \
    	  "/projects/pfenninggroup/machineLearningForComputationalBiology/halLiftover_chains/data/raw_data/fasta/Mus_musculus.fa" \
    	  -out_file "/home/azstephe/regression_liver/data/model_outputs/${modelname}/outputs.csv" \
    	  --write_csv
"
