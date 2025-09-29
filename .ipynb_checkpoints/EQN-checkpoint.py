#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Processes BED files using Extended Quantile Normalization (EQN) Algorithm.

This script takes a BED file name and a base directory, then processes
several versions of the file (from 'log', 'log_sorted', 'quantile_norm'
subdirectories).

It splits the peaks into an "included" and "excluded" set. For each peak
in the excluded set, it finds the closest peak value from the included set and
assigns its quantile-normalized signal value. The final, combined data
is saved to a 'eqn' subdirectory.
"""

import pandas as pd
import argparse
import os

def find_closest(row, sorted_df):
    """
    Finds the name of the peak in sorted_df with the value in column 4
    that is closest to the value in the given row's column 4.
    """
    # Get idx with minimum difference in signal
    closest_idx = (sorted_df[4] - row[4]).abs().idxmin()
    return sorted_df.loc[closest_idx, 3]

def process_data(filename, base_dir):
    """Main processing logic."""
    print(f"Processing file: {filename} in directory: {base_dir}")

    # Load data - assumed log, log_sorted, and quantile normalized data already processed and created
    try:
        log_all = pd.read_csv(os.path.join(base_dir, 'log', filename), header=None, sep='\t')
        log_sorted = pd.read_csv(os.path.join(base_dir, 'log_sorted', filename), header=None, sep='\t')
        qn_included = pd.read_csv(os.path.join(base_dir, 'quantile_norm', filename), header=None, sep='\t')
    except FileNotFoundError as e:
        print(f"Error loading files: {e}")
        print("Please ensure the directory structure ('log', 'log_sorted', 'quantile_norm') exists within the base directory.")
        return

    # Split data into included/exluded sets
    # Hardcoded based on smallest sample (e.g. pig)
    split_index = 20615
    excluded = log_all[split_index:].copy()
    
    # Create the included set by filtering the sorted log file
    sorted_included = log_sorted[~log_sorted[3].isin(excluded[3])].reset_index(drop=True)

    # Find closest peaks for the excluded set
    print("Finding closest peaks for the excluded set... (this may take a while)")
    excluded['closest_col3'] = excluded.apply(lambda row: find_closest(row, sorted_included), axis=1)
    print("Done finding closest peaks.")

    # Merge to get the quantile normalized signal for the excluded peaks
    qn_excluded = excluded.merge(qn[[3, 4]], left_on='closest_col3', right_on=3, how='left')

    # Clean up columns after the merge
    qn_excluded.drop(columns=['closest_col3', '3_y', '4_x'], inplace=True)
    qn_excluded.rename(columns={'3_x': 3, '4_y': 4}, inplace=True)

    # Concatenate qn values from included and excluded sets
    eqn = pd.concat([qn_included, qn_excluded]).reset_index(drop=True)

    output_dir = os.path.join(base_dir, 'eqn')
    os.makedirs(output_dir, exist_ok=True) # Create directory if it doesn't exist
    output_path = os.path.join(output_dir, filename)
    
    eqn.to_csv(output_path, header=None, index=False, sep='\t')
    print(f"\nSuccessfully saved final data to: {output_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process BED files to create an extended quantile normalized dataset.")
    
    parser.add_argument('filename', type=str, 
                        help="The name of the BED file to process (e.g., 'macaque_liver_pos_ALL.bed').")
    
    parser.add_argument('directory', type=str, 
                        help="The base directory containing the 'log', 'log_sorted', and 'quantile_norm' subdirectories.")
    
    args = parser.parse_args()
    
    process_data(args.filename, args.directory)