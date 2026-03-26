from typing import Any
import numpy as np
import tensorflow as tf
import tensorflow_hub as hub
from Bio import SeqIO
import os
from pathlib import Path
import json

def one_hot_encode(sequence: str,
                   alphabet: str = 'ACGT',
                   neutral_alphabet: str = 'N',
                   neutral_value: Any = 0,
                   dtype=np.float32) -> np.ndarray:
  """One-hot encode sequence."""
  def to_uint8(string):
    return np.frombuffer(string.encode('ascii'), dtype=np.uint8)
  hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
  hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
  hash_table[to_uint8(neutral_alphabet)] = neutral_value
  hash_table = hash_table.astype(dtype)
  return hash_table[to_uint8(sequence)]

def prepare_input(sequence: str):
    target_length = 393216
    if len(sequence) != target_length:
        print("input length != 393216")
        return

    # Use your encoding function
    encoded = one_hot_encode(sequence)
    
    # Add the Batch Dimension: (393216, 4) -> (1, 393216, 4)
    return encoded[np.newaxis, ...]
    
def process_fasta_to_batch(fasta_path):
    encoded_list = []
    
    # 1. Iterate through the FASTA file
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        
        # 2. Use your prepare_input logic
        # Note: prepare_input already adds the (1, 393216, 4) dimension
        encoded_seq = prepare_input(seq)
        
        if encoded_seq is not None:
            encoded_list.append(encoded_seq)
    
    # 3. Stack them into a single batch
    # Since prepare_input returns (1, L, 4), we use concatenate on axis 0
    if encoded_list:
        batch_array = np.concatenate(encoded_list, axis=0)
        return batch_array
    else:
        return None

enformer_model = hub.load("https://www.kaggle.com/models/deepmind/enformer/TensorFlow2/enformer/1").model
ENFORMER_INPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/enformer_inputs/")
ENFORMER_OUTPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/enformer_outputs/")
target_track = 205

dir_list = ['log_pos', 'log_test3', 'neg']
# 'log_test1', 'log_test2'
# for subdir_name in os.listdir(ENFORMER_INPUTS):

for subdir_name in dir_list:
    subdir_path = ENFORMER_INPUTS / subdir_name
    if not subdir_path.is_dir():
        continue

    for fasta_path in subdir_path.glob("*.fa"):
        print(f"--- Processing: {fasta_path} ---", flush=True)

        # 1. Setup paths
        stem = fasta_path.stem.replace("_500bp", "")
        out_subdir = ENFORMER_OUTPUTS / subdir_name
        out_subdir.mkdir(parents=True, exist_ok=True)
        
        save_path = out_subdir / f"{stem}_enformer_{target_track}_preds.npy"

        # 2. Skip check: check the OUTPUT file, not the input
        if save_path.exists() and save_path.stat().st_size > 0:
            print(f"Skipping: {save_path.name} already exists.", flush=True)
            continue

#for subdir_name in dir_list:
#    subdir_path = ENFORMER_INPUTS / subdir_name
#    if not subdir_path.is_dir(): continue
#
#    for fasta_path in subdir_path.glob("*.fa"):
#        print(f"--- Processing: {fasta_path} ---", flush=True)
#
#	stem = fasta_path.stem.replace("_500bp", "")
#        out_subdir = ENFORMER_OUTPUTS / subdir_name
#        out_subdir.mkdir(parents=True, exist_ok=True)
#
#        save_path = out_subdir / f"{stem}_enformer_{target_track}_preds.npy#"
#
#	if save_path.exists() and save_path.stat().st_size > 0:
#	    print(f"Skipping: {fasta_path.name} already exists and is non-empty.", flush=True)
#            continue

        results = [] 
        x = process_fasta_to_batch(str(fasta_path)) 
        
        for i in range(len(x)):
            # Grab one sequence: Shape (1, 196608, 4)
            single_input = x[i : i+1]
            
            # Predict
            pred = enformer_model.predict_on_batch(single_input)
            
            track_subset = pred['mouse'][0, :, target_track]
            
            # Convert to numpy immediately to save GPU memory
            results.append(track_subset.numpy() if hasattr(track_subset, 'numpy') else track_subset)
            
            if i % 10 == 0:
                print(f"  Sequence {i}/{len(x)}", flush=True)

        # 3. CONVERT TO ARRAY: Shape (Num_Peaks, 896)
        final_output = np.array(results)
        
        #stem = fasta_path.stem.replace("_500bp", "")
        #out_subdir = ENFORMER_OUTPUTS / subdir_name
        #out_subdir.mkdir(parents=True, exist_ok=True)
        
        #save_path = out_subdir / f"{stem}_enformer_{target_track}_preds.npy"
        np.save(save_path, final_output)

