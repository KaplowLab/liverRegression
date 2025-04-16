import pandas as pd
import numpy as np
import sys

def log_transform_bed(file_path, output_path):
    # Load the BED file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None)

    # Ensure the file has at least 7 columns (BED6 + signal column)
    if df.shape[1] < 5:
        raise ValueError("The input file must have at least 7 columns.")

    # Extract the signal column (7th column, index 6)
    signals = df[4].values

    # Apply log transformation (adding a small constant to avoid log(0))
    df[4] = np.log1p(signals)

    # Save the modified DataFrame to a new BED file
    df.to_csv(output_path, sep='\t', header=False, index=False)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = '/'.join(input_file.split('/')[:-1]) + "/log/" + input_file.split('/')[-1]
    log_transform_bed(input_file, output_file)
