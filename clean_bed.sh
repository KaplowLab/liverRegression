CRIPT CONFIGURATION
# ==============================================================================
# Please define the paths to your files here for easy editing.

# --- INPUT FILES ---
BED_FILE_ORIGINAL="/home/azstephe/liverRegression/regression_liver/data/test_splits_2kb/log_pos_LiuAll/rat_liver_TEST_500bp.bed"
GENOME_FILE="/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatGenome/rn6.fa"
# --- OUTPUT FILES ---
# The script will create this final, cleaned BED file.
CLEAN_BED_FILE="/home/azstephe/liverRegression/regression_liver/data/test_splits_2kb/log_pos_LiuAll/rat_liver_TEST_500bp_clean.bed"

# You can leave these temporary filenames as they are.
TEMP_FASTA_OUTPUT="temp_fasta_output.fa"
TEMP_WARNING_LOG="temp_warnings.log"
TEMP_COORDS_TO_REMOVE="temp_coords_to_remove.bed"


# ==============================================================================
# MAIN SCRIPT LOGIC
# ==============================================================================
# Exit immediately if any command fails to prevent corrupted results.
set -e

echo "--- Step 1 of 4: Running bedtools and capturing warnings... ---"
# Run the getfasta command. We redirect its standard error stream (2>), which
# contains the warnings, into our temporary log file.
bedtools getfasta -fo "$TEMP_FASTA_OUTPUT" -fi "$GENOME_FILE" -bed "$BED_FILE_ORIGINAL" 2> "$TEMP_WARNING_LOG"
echo "Warnings have been successfully captured in '$TEMP_WARNING_LOG'"


echo "\n--- Step 2 of 4: Extracting problematic coordinates from warnings... ---"
# Parse the log file to create a BED file of the coordinates that caused errors.
# The -F flag sets the field separators to any of '(', ')', ':', or '-'.
awk -F'[():-]' '{print $2 "\t" $3 "\t" $4}' "$TEMP_WARNING_LOG" > "$TEMP_COORDS_TO_REMOVE"
echo "Problematic coordinates have been extracted to '$TEMP_COORDS_TO_REMOVE'"


echo "\n--- Step 3 of 4: Filtering the original BED file... ---"
# Use bedtools subtract to remove the problematic regions from the original BED file.
# The '-a' file is the source, and we subtract the regions from the '-b' file.
bedtools subtract -a "$BED_FILE_ORIGINAL" -b "$TEMP_COORDS_TO_REMOVE" > "$CLEAN_BED_FILE"
echo "Original BED file has been cleaned."


echo "\n--- Step 4 of 4: Cleaning up temporary files... ---"
# Remove the intermediate files that are no longer needed.
#rm "$TEMP_WARNING_LOG" "$TEMP_COORDS_TO_REMOVE" "$TEMP_FASTA_OUTPUT"
echo "Temporary files have been removed."


echo "\n--------------------------------------------------"
echo "Success! Your clean BED file is ready at:"
echo "$CLEAN_BED_FILE"
echo "--------------------------------------------------"
