#!/bin/bash

# --- Configuration ---
# Set the input and output directories.
# Make sure to change these to your actual directory paths.
INPUT_DIR=/home/azstephe/liverRegression/regression_liver/data/test_splits/$1
OUTPUT_DIR=/home/azstephe/liverRegression/regression_liver/data/test_splits_2kb/$1

# --- Main Script ---

# Exit immediately if a command fails.
set -e

# Create the output directory if it doesn't already exist.
# The '-p' flag prevents an error if the directory is already there.
echo "Ensuring output directory '$OUTPUT_DIR' exists..."
mkdir -p "$OUTPUT_DIR"

# Check if the input directory is empty
if [ -z "$(ls -A "$INPUT_DIR")" ]; then
   echo "Warning: Input directory '$INPUT_DIR' is empty. No files to process."
   exit 0
fi


# Loop through every file in the input directory.
# Using quotes around "$INPUT_DIR/*" handles filenames with spaces.
for input_file in "$INPUT_DIR"/*
do
  # Get just the filename from the full path (e.g., 'path/to/file.txt' -> 'file.txt')
  filename=$(basename "$input_file")

  # Construct the full path for the output file.
  output_file="$OUTPUT_DIR/$filename"

  # Print a status message to the console.
  echo "Processing '$input_file' -> '$output_file'"

  # Run your awk command on the input file and redirect the output to the new file
  awk 'BEGIN{OFS="\t"} {print $1, $2 - 750, $3 + 750, $4, $5}' "$input_file" | awk '$2 >= 0' > "$output_file"
done

echo "---------------------"
echo "All files processed successfully to $OUTPUT_DIR"
