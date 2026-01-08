#!/bin/bash
#
# Usage: bash setup.sh <cluster_name>, where <cluster_name> is "lane" or "bridges"

# Check arguments

if [ "$CONDA_DEFAULT_ENV" != "keras2-tf27" ]; then
        source activate keras2-tf27
fi

MODEL=$1
VAL1=$2
VAL2=$3
VAL3=$4
GENOME=$5
OUTPUT=$6
OUTPUT_NAME=$7
python scripts/get_activations.py -model $MODEL -in_files $VAL1 $VAL2 $VAL3 -in_genomes $GENOME $GENOME $GENOME -out_file $OUTPUT/$OUTPUT_NAME --write_csv


