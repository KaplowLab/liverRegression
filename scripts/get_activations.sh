#!/bin/bash

# Change the -p [partition] in submit_job() helper function

# Usage
if [ "$#" -ne 2 ]; then
    echo "Usage: bash get_activations.sh <PROJECT_DIR> <RUN_DIR>"
    echo "Example: bash get_activations.sh /home/user/project /path/to/log/run-20250225-bdbi7l3n"
    exit 1
fi

PROJECT_DIR=$1
RUN_DIR=$2

# 2. Parsing Run Info
# Extracts 'run-20250225-bdbi7l3n' from the path
RUN_NAME=$(basename "$RUN_DIR")
# Extracts 'bdbi7l3n' (everything after the last hyphen)
ID=${RUN_NAME##*-}

# Construct internal paths
MODEL_PATH="$RUN_DIR/files/model-best.h5"
OUTPUT_BASE="${PROJECT_DIR}/data/model_outputs/${ID}"

# Create output directory
mkdir -p "$OUTPUT_BASE"

# Environment Setup
if [ "${CONDA_DEFAULT_ENV:-}" != "keras2-tf27" ]; then
    eval "$(conda shell.bash hook)"
    conda activate keras2-tf27
fi

# HELPER FUNCTION
# Usage: submit_job <2|3> <species_label> <output_filename> <genome_path> <input1> <input2> [input3]
submit_job() {
    local mode=$1; local label=$2; local out_name=$3; local genome=$4; local v1=$5; local v2=$6; local v3=${7:-""}
    local script="get_activations_helper${mode}.sh"
    local job_name="${label}_${OUTPUT_BASE: -3}"
    local log_file="${label}_${OUTPUT_BASE: -3}.o"

    echo "Submitting $label ($mode-way)..."

    if [ "$mode" -eq 2 ]; then
        sbatch --mem=4G -o "$log_file" -J "$job_name" -p gpu -n 1 --gres gpu:1 "$script" "$MODEL" "$v1" "$v2" "$genome" "$OUTPUT_BASE" "$out_name"
    else
        sbatch --mem=4G -o "$log_file" -J "$job_name" -p gpu -n 1 --gres gpu:1 "$script" "$MODEL" "$v1" "$v2" "$v3" "$genome" "$OUTPUT_BASE" "$out_name"
    fi
}

# Update with user genomes paths (can be downloaded from UCSC)
MM10="${PROJECT_DIR}/data/genomes/mm10.fa"
RHEMAC8="${PROJECT_DIR}/data/genomes/rheMac8.fa"
RN6="${PROJECT_DIR}/data/genomes/rn6.fa"
BTAU="${PROJECT_DIR}/data/genomes/Btau.fna"
SUS_SCR="${PROJECT_DIR}/data/genomes/susScr3.fa"

# ==============================================================================
# MOUSE
# ==============================================================================
submit_job 2 "moTR" "activations_mouse_TRAIN.csv" "$MM10" \
    "${PROJECT_DIR}/data/train_splits/log_pos/mouse_liver_TRAIN_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/neg/mouse_liver_TRAIN_500bp.bed"

submit_job 2 "moV" "activations_mouse_VAL.csv" "$MM10" \
    "${PROJECT_DIR}/data/val_splits/log_pos/mouse_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/neg/mouse_liver_VAL_500bp.bed"

submit_job 2 "moTE" "activations_mouse_TEST.csv" "$MM10" \
    "${PROJECT_DIR}/data/test_splits/log_pos/mouse_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/neg/mouse_liver_TEST_500bp.bed"

# ==============================================================================
# MACAQUE
# ==============================================================================
submit_job 2 "macTR" "activations_macaque_TRAIN.csv" "$RHEMAC8" \
    "${PROJECT_DIR}/data/train_splits/log_pos/macaque_liver_TRAIN_500bp.bed" \
    "${PROJECT_DIR}/data/train_splits/neg/macaque_liver_TRAIN_500bp.bed"

submit_job 3 "macV" "activations_macaque_VAL.csv" "$RHEMAC8" \
    "${PROJECT_DIR}/data/val_splits/val1/macaque_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val2/macaque_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val3/macaque_liver_VAL_500bp.bed"

submit_job 3 "macTE" "activations_macaque_TEST.csv" "$RHEMAC8" \
    "${PROJECT_DIR}/data/test_splits/test1/macaque_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test2/macaque_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test3/macaque_liver_TEST_500bp.bed"

submit_job 2 "macVpn" "activations_macaque_VAL_pos_neg.csv" "$RHEMAC8" \
    "${PROJECT_DIR}/data/val_splits/log_pos/macaque_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/neg/macaque_liver_VAL_500bp.bed"

submit_job 2 "macTEpn" "activations_macaque_TEST_pos_neg.csv" "$RHEMAC8" \
    "${PROJECT_DIR}/data/test_splits/log_pos/macaque_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/neg/macaque_liver_TEST_500bp.bed"

# ==============================================================================
# RAT (RN6)
# ==============================================================================
submit_job 2 "ratTR" "activations_rat_TRAIN.csv" "$RN6" \
    "${PROJECT_DIR}/data/train_splits/log_pos/rat_liver_TRAIN_500bp.bed" \
    "${PROJECT_DIR}/data/train_splits/neg/rat_liver_TRAIN_500bp.bed"

submit_job 3 "ratV" "activations_rat_VAL.csv" "$RN6" \
    "${PROJECT_DIR}/data/val_splits/val1/rat_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val2/rat_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val3/rat_liver_VAL_500bp.bed"

submit_job 3 "ratTE" "activations_rat_TEST.csv" "$RN6" \
    "${PROJECT_DIR}/data/test_splits/test1/rat_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test2/rat_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test3/rat_liver_TEST_500bp.bed"

submit_job 2 "ratVpn" "activations_rat_VAL_pos_neg.csv" "$RN6" \
    "${PROJECT_DIR}/data/val_splits/log_pos/rat_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/neg/rat_liver_VAL_500bp.bed"

submit_job 2 "ratTEpn" "activations_rat_TEST_pos_neg.csv" "$RN6" \
    "${PROJECT_DIR}/data/test_splits/log_pos/rat_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/neg/rat_liver_TEST_500bp.bed"

# ==============================================================================
# COW (BTAU)
# ==============================================================================
submit_job 2 "cowTR" "activations_cow_TRAIN.csv" "$BTAU" \
    "${PROJECT_DIR}/data/train_splits/log_pos/cow_liver_TRAIN_500bp.bed" \
    "${PROJECT_DIR}/data/train_splits/neg/cow_liver_TRAIN_500bp.bed"

submit_job 3 "cowV" "activations_cow_VAL.csv" "$BTAU" \
    "${PROJECT_DIR}/data/val_splits/val1/cow_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val2/cow_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val3/cow_liver_VAL_500bp.bed"

submit_job 3 "cowTE" "activations_cow_TEST.csv" "$BTAU" \
    "${PROJECT_DIR}/data/test_splits/test1/cow_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test2/cow_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test3/cow_liver_TEST_500bp.bed"

submit_job 2 "cowVpn" "activations_cow_VAL_pos_neg.csv" "$BTAU" \
    "${PROJECT_DIR}/data/val_splits/log_pos/cow_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/neg/cow_liver_VAL_500bp.bed"

submit_job 2 "cowTEpn" "activations_cow_TEST_pos_neg.csv" "$BTAU" \
    "${PROJECT_DIR}/data/test_splits/log_pos/cow_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/neg/cow_liver_TEST_500bp.bed"

# ==============================================================================
# PIG (SUS_SCR)
# ==============================================================================
submit_job 2 "pigTR" "activations_pig_TRAIN.csv" "$SUS_SCR" \
    "${PROJECT_DIR}/data/train_splits/log_pos/pig_liver_TRAIN_500bp.bed" \
    "${PROJECT_DIR}/data/train_splits/neg/pig_liver_TRAIN_500bp.bed"

submit_job 3 "pigV" "activations_pig_VAL.csv" "$SUS_SCR" \
    "${PROJECT_DIR}/data/val_splits/val1/pig_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val2/pig_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/log_val3/pig_liver_VAL_500bp.bed"

submit_job 3 "pigTE" "activations_pig_TEST.csv" "$SUS_SCR" \
    "${PROJECT_DIR}/data/test_splits/test1/pig_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test2/pig_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/log_test3/pig_liver_TEST_500bp.bed"

submit_job 2 "pigVpn" "activations_pig_VAL_pos_neg.csv" "$SUS_SCR" \
    "${PROJECT_DIR}/data/val_splits/log_pos/pig_liver_VAL_500bp.bed" \
    "${PROJECT_DIR}/data/val_splits/neg/pig_liver_VAL_500bp.bed"

submit_job 2 "pigTEpn" "activations_pig_TEST_pos_neg.csv" "$SUS_SCR" \
    "${PROJECT_DIR}/data/test_splits/log_pos/pig_liver_TEST_500bp.bed" \
    "${PROJECT_DIR}/data/test_splits/neg/pig_liver_TEST_500bp.bed"
