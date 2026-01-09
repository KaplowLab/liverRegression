import pandas as pd
import numpy as np
import sys
import os

def log_transform_bed(file_path, output_path):
    """
    Reads a narrowPeak file, applies a log(x+1) transformation to the 
    signal column, and saves output to a new directory.
    
    Parameters:
    file_path (str): Path to the input .narrowPeak file.
    output_path (str): Path where the transformed file will be saved.
    """
    # Load the BED file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None)

    if df.shape[1] < 7:
        raise ValueError("The input file must have at least 7 columns")

    signals = df[6].values

    # Apply log1p transformation: log(1 + x)
    df[6] = np.log1p(signals)

    df.to_csv(output_path, sep='\t', header=False, index=False)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file>")
    else:

    	input_file = sys.argv[1]
    
    	# Construct output path: places the file in a /log/ subfolder within the same directory
    	# Example: data/peaks.bed -> data/log/peaks.bed
    	dir_name = os.path.dirname(input_file)
    	file_name = os.path.basename(input_file)
    	log_dir = os.path.join(dir_name, "log_transformed")

    	os.makedirs(log_dir, exist_ok=True)
    	output_file = os.path.join(log_dir, file_name)

    	# Run the transformation
    	log_transform_bed(input_file, output_file)
