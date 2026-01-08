#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Quantile normalizes input files

This script takes multiple space separated bed or narrowPeak file name paths. Assumes the files are the same length. Change the value of i in the load_bed_file function to denote the column containing the signal value (5 for BED file; 7 for narrowPeak file).

Saves quantile normalized output files in "quantile_norm" subdirectory

Quantile normalizes with no ties for rank.
"""

import numpy as np
from scipy.stats import rankdata
import pandas as pd
import sys
from pathlib import Path
import os

def rank_column_with_ties(col):
    """Custom function for ranking that assigns different ranks to tied values"""
    sorted_indices = np.argsort(col)
    ranks = np.zeros_like(col, dtype=int)
    
    # Assign ranks based on sorted order without considering ties
    for i in range(len(col)):
        ranks[sorted_indices[i]] = i + 1 
    return ranks

def quantile_normalize(values):
    """Quantile normalizes 2D array"""
    sorted_values = np.sort(values, axis=1)
    ranked = np.apply_along_axis(rank_column_with_ties, 1, values)

    # Calculate mean value at each rank
    quantiles = np.mean(sorted_values, axis=0)

    # Vectorize replacement function
    replacement_map = {i+1: quantiles[i] for i in range(len(quantiles))}
    vectorized_replace = np.vectorize(lambda x: replacement_map.get(x, x))

    # Replace ranks with corresponding mean values
    normalized_values = vectorized_replace(ranked)
    return normalized_values

def load_bed_file(input_file, i):
    df = pd.read_csv(input_file, sep='\t', header=None)
    if df.shape[1] < i:
        raise ValueError(f'The input file must have at least {i} columns.')
    return df

def save_bed_file(df, output_file):
    df.to_csv(output_file, sep='\t', header=False, index=False)

def main(input_files, output_files, i=5):
    """Saves quantile normalized data to output_files"""
    dfs = [load_bed_file(file, i) for file in input_files]
    
    # 2D array - each df occupies a row
    all_values = np.vstack([df[i-1].values for df in dfs])
    
    # Quantile normalize the ith column
    all_values_qn = quantile_normalize(all_values)
    
    # Return the normalized values back into individual dfs
    index = 0
    for i, df in enumerate(dfs):
        size = df.shape[0]
        df[i-1] = all_values_qn[i]
        save_bed_file(df, output_files[i])

def new_path_dir(file_paths):
    """Creates path name for outfiles"""
    updated_paths = []
    for file_path in file_paths:
        path_obj = Path(file_path)

        parent_dir = path_obj.parent
        
        new_dir = Path(*parent_dir.parts, 'quantile_norm')
        updated_path = new_dir / path_obj.name
        
        new_dir.mkdir(parents=True, exist_ok=True)
        
        updated_paths.append(str(updated_path))
    return updated_paths
    
if __name__ == "__main__":
    """
    Usage: Follow python command with input files separated by a space
    Example: python quantile_normalize.py macaque.bed rat.bed mouse.bed
    """
    infiles = sys.argv[1:]

    outfiles = new_path_dir(infiles)
    
    main(infiles, outfiles)
    

