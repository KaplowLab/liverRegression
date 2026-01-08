#!/bin/bash

# Description: Gets model activations for all train, validation, and test evaluation sets.
# NOTE: $MODEL path is hardcoded so feel free to change.
# Usage: bash get_activations.sh <model dir>
# Example: bash get_activations.sh run-20250225_222915-bdbi7l3n

if [ "$CONDA_DEFAULT_ENV" != "keras2-tf27" ]; then
        source activate keras2-tf27
fi

RUN=$1
ID=${RUN##*-}
MODEL=/home/azstephe/liverRegression/regression_liver/data/best_models/log/$RUN/files/model-best.h5
OUTPUT=/home/azstephe/liverRegression/regression_liver/data/model_outputs/${ID}_FINAL

if [ ! -d $OUTPUT ]; then
   mkdir $OUTPUT
fi

# Mouse posTrain, posVal, and negVal activations

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/mouse_liver_TRAIN_500bp.narrowPeak

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/mouse_liver_VAL_500bp.narrowPeak

VAL3=/home/azstephe/regression_liver/data/splits/negatives/nonMouse_liver_andRat_andCow_andPig_andMacaque_VAL_500bp.bed

GENOME=/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/mm10.fa

OUTPUT_NAME=activations_mouse_TRAIN_VAL.csv

sbatch --mem=4G -o mo_activ.o -J mo3_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos/mouse_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/neg/mouse_liver_TEST_500bp.bed

GENOME=/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/mm10.fa 

OUTPUT_NAME=activations_mouse_TEST.csv 

sbatch --mem=4G -o mo_activ.o -J mo_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

# MACAQUE 

# Pos and Neg Train
VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/macaque_liver_TRAINONLY.narrowPeak 

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonMacaque_liver_andRat_andCow_andPig_TRAIN_500bp.bed 

GENOME=/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueGenome/rheMac8.fa

OUTPUT_NAME=activations_macaque_TRAIN.csv

sbatch --mem=4G -o mac_activ.o -J mac_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

# Val 1,2,3
VAL1=/home/azstephe/liverRegression/regression_liver/data/val_splits/val1/macaque_liver_VAL_500bp.bed 

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/log_val2/macaque_liver_VAL.narrowPeak 

VAL3=/home/azstephe/liverRegression/regression_liver/data/splits/log_val3/macaque_liver_VAL.narrowPeak 

OUTPUT_NAME=activations_macaque_VAL.csv

sbatch --mem=4G -o mac_activ.o -J mac_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME

# Pos and Neg Val
VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/macaque_liver_VAL.narrowPeak

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonMacaque_liver_andRat_andCow_andPig_VAL_500bp.bed 

OUTPUT_NAME=activations_macaque_VAL_orthologs.csv

sbatch --mem=4G -o mac_activ.o -J mac_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

# Pos and Neg Test
VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/macaque_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/neg/macaque_liver_TEST_500bp.bed

OUTPUT_NAME=activations_macaque_TEST_orthologs.csv

sbatch --mem=4G -o mac_activ.o -J mac_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

# Test 1,2,3
VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/macaque_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/macaque_liver_TEST_500bp.bed

VAL3=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/macaque_liver_TEST_500bp.bed

GENOME=/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueGenome/rheMac8.fa

OUTPUT_NAME=activations_macaque_TEST.csv

sbatch --mem=4G -o mac_activ.o -J mac_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME

# RAT

# Pos and Neg Train
VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/rat_liver_TRAINONLY.narrowPeak

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonRat_liver_andMacaque_andCow_andPig_TRAIN_500bp.bed

GENOME=/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatGenome/rn6.fa

OUTPUT_NAME=activations_rat_TRAIN.csv

sbatch --mem=4G -o rat_activ.o -J rat_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/rat_liver_VAL.narrowPeak 

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonRat_liver_andMacaque_andCow_andPig_VAL_500bp.bed 

OUTPUT_NAME=activations_rat_VAL_orthologs.csv

sbatch --mem=4G -o rat_activ.o -J rat_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/val_splits/val1/rat_liver_VAL_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/log_val2/rat_liver_VAL.narrowPeak

VAL3=/home/azstephe/liverRegression/regression_liver/data/splits/log_val3/rat_liver_VAL.narrowPeak

OUTPUT_NAME=activations_rat_VAL.csv

sbatch --mem=4G -o rat_activ.o -J rat_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/rat_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/neg/rat_liver_TEST_500bp.bed

OUTPUT_NAME=activations_rat_TEST_orthologs.csv  

sbatch --mem=4G -o rat_activ.o -J rat_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME
#sbatch --mem=4G -o ratb_activ.o -J ratb_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2_baye.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT_BAYE $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/rat_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/rat_liver_TEST_500bp.bed

VAL3=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/rat_liver_TEST_500bp.bed

OUTPUT_NAME=activations_rat_TEST.csv

sbatch --mem=4G -o rat_activ.o -J rat_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME
#sbatch --mem=4G -o ratb_activ.o -J ratb_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3_baye.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT_BAYE $OUTPUT_NAME

# cow

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/cow_liver_TRAINONLY.narrowPeak

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonCow_liver_andMacaque_andRat_andPig_TRAIN_500bp.bed

GENOME=/home/azstephe/regression_liver/data/splits/cowMouse/GCF_000003205.7_Btau_5.0.1_genomic_chromName.fna

OUTPUT_NAME=activations_cow_TRAIN.csv

sbatch --mem=4G -o cow_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/cow_liver_VAL.narrowPeak 

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonCow_liver_andMacaque_andRat_andPig_VAL_500bp.bed 

OUTPUT_NAME=activations_cow_VAL_orthologs.csv

sbatch --mem=4G -o cow_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/val_splits/val1/cow_liver_VAL_500bp.bed 

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/log_val2/cow_liver_VAL.narrowPeak 

VAL3=/home/azstephe/liverRegression/regression_liver/data/splits/log_val3/cow_liver_VAL.narrowPeak

OUTPUT_NAME=activations_cow_VAL.csv

sbatch --mem=4G -o cow_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/cow_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/neg/cow_liver_TEST_500bp.bed

OUTPUT_NAME=activations_cow_TEST_orthologs.csv  

sbatch --mem=4G -o cow_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME
#sbatch --mem=4G -o cowb_activ.o -J cowb_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2_baye.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT_BAYE $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/cow_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/cow_liver_TEST_500bp.bed

VAL3=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/cow_liver_TEST_500bp.bed

OUTPUT_NAME=activations_cow_TEST.csv

sbatch --mem=4G -o cow_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME
#sbatch --mem=4G -o cowb_activ.o -J cowb_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3_baye.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT_BAYE $OUTPUT_NAME

# pig 

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/pig_liver_TRAINONLY.narrowPeak

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonPig_liver_andMacaque_andRat_andCow_TRAIN_500bp.bed 

GENOME=/data/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/susScr3.fa

OUTPUT_NAME=activations_pig_TRAIN.csv

sbatch --mem=16G -o pig_activ.o -J pig_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/logPos/pig_liver_VAL.narrowPeak 

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/negatives/nonPig_liver_andMacaque_andRat_andCow_VAL_500bp.bed 

OUTPUT_NAME=activations_pig_VAL_orthologs.csv

sbatch --mem=4G -o pig_activ.o -J pig_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/val_splits/val1/pig_liver_VAL_500bp.bed 

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/log_val2/pig_liver_VAL.narrowPeak 

VAL3=/home/azstephe/liverRegression/regression_liver/data/splits/log_val3/pig_liver_VAL.narrowPeak 

OUTPUT_NAME=activations_pig_VAL.csv

sbatch --mem=4G -o pig_activ.o -J pig_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/pig_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/neg/pig_liver_TEST_500bp.bed

OUTPUT_NAME=activations_pig_TEST_orthologs.csv

sbatch --mem=4G -o pig_activ.o -J pig_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/pig_liver_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/pig_liver_TEST_500bp.bed

VAL3=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/pig_liver_TEST_500bp.bed

GENOME=/data/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/susScr3.fa

OUTPUT_NAME=activations_pig_TEST.csv

sbatch --mem=4G -o pig_activ.o -J pig_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ3.sh $MODEL $VAL1 $VAL2 $VAL3 $GENOME $OUTPUT $OUTPUT_NAME

# pig + cow

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_TRAIN_chromName_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_TRAIN_chromName_500bp.bed

GENOME=/home/azstephe/regression_liver/data/splits/cowMouse/GCF_000003205.7_Btau_5.0.1_genomic_chromName.fna

OUTPUT_NAME=activations_cow_pig_TRAIN.csv

sbatch --mem=4G -o cow_pig_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_VAL_chromName_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_VAL_chromName_500bp.bed

OUTPUT_NAME=activations_cow_pig_VAL.csv

sbatch --mem=4G -o cow_pig_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

VAL1=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test4/cow_pig_liver_pos_TEST_500bp.bed

VAL2=/home/azstephe/liverRegression/regression_liver/data/test_splits/log_test5/cow_pig_liver_neg_TEST_500bp.bed

OUTPUT_NAME=activations_cow_pig_TEST.csv

sbatch --mem=4G -o cow_pig_activ.o -J cow_"${OUTPUT: -3}" -p gpu -n 1 --gres gpu:1 get_activ2.sh $MODEL $VAL1 $VAL2 $GENOME $OUTPUT $OUTPUT_NAME

