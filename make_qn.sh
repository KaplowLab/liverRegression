#!/bin/bash

# This script automates the process of running filterPeakName.py
# for a predefined set of directories and species.

# --- Configuration ---
# An array to hold the directory names
DIRS=("log_pos" "log_pos_LiuAll" "log_test2" "log_test3" "log_test4")

# An array to hold the species names
SPECIES_LIST=("mouse" "macaque" "rat" "cow" "pig")

# --- Main Execution ---
# Loop over each directory in the DIRS array
for DIR in "${DIRS[@]}"; do
    echo "--- Processing directory: ${DIR} ---"

    # Before processing species, ensure the output directory exists
    # The 'mkdir -p' command creates the directory and any parent directories if needed
    OUTPUT_PARENT_DIR="/home/azstephe/liverRegression/regression_liver/data/test_splits_qn/${DIR}"
    mkdir -p "${OUTPUT_PARENT_DIR}"

    # Loop over each species in the SPECIES_LIST array
    for SPECIES in "${SPECIES_LIST[@]}"; do
        echo "  -> Running for species: ${SPECIES}"

        # Define the input and output file paths using the current DIR and SPECIES
        UNFILTERED_PEAK_FILE="/home/azstephe/liverRegression/regression_liver/data/sorted_log_20615/quantile_norm/${SPECIES}_liver_pos_ALL.bed"
        PEAK_LIST_FILE="/home/azstephe/liverRegression/regression_liver/data/test_splits_2kb/${DIR}/${SPECIES}_liver_TEST_500bp.bed"
        OUTPUT_FILE="${OUTPUT_PARENT_DIR}/${SPECIES}_liver_TEST_500bp.bed"

        # Run the Python script with the specified parameters
        # The backslashes (\) allow us to break the command into multiple lines for readability
        python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py \
            --unfilteredPeakFileName "${UNFILTERED_PEAK_FILE}" \
            --peakListFileName "${PEAK_LIST_FILE}" \
            --peakNameCol 3 \
            --outputFileName "${OUTPUT_FILE}"

    done
    echo "--- Finished processing directory: ${DIR} ---"
    echo "" # Add a newline for cleaner output
done

echo "All jobs completed! 🎉"
