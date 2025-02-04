import numpy as np
import pandas as pd
from scipy.stats import rankdata
import sys

def quantile_normalize(values):
    """Quantile normalization of a vector."""
    sorted_idx = np.argsort(values)
    sorted_values = np.sort(values)
    quantiles = np.mean(sorted_values)
    return np.array([quantiles for i in range(len(values))])[sorted_idx.argsort()]

def main(input_file, output_file):
    # Load the BED file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None)

    # Ensure the file has at least 7 columns (BED6 + signal column)
    if df.shape[1] < 7:
        raise ValueError("The input file must have at least 7 columns.")

    # Log-transform the values in the seventh column (log base 10)
    df[6] = np.log10(df[6].replace(0, np.nan))  # Avoid log(0) issues by replacing 0 with NaN

    # Drop NaN values (optional, if there are zeros in your data)
    df = df.dropna(subset=[6])

    # Quantile normalize the log-transformed values in the seventh column
    df[6] = quantile_normalize(df[6].values)

    # Save the result to a new BED file
    df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"Log-transformation and quantile normalization completed. Output saved to {output_file}.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Run main function
    main(input_file, output_file)

