from typing import Any
import numpy as np
import tensorflow as tf
import tensorflow_hub as hub
from Bio import SeqIO
import os
from pathlib import Path
import json

# from https://github.com/google-deepmind/deepmind-research/tree/master/enformer
def one_hot_encode(sequence: str,
                   alphabet: str = 'ACGT',
                   neutral_alphabet: str = 'N',
                   neutral_value: Any = 0,
                   dtype=np.float32) -> np.ndarray:
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

    encoded = one_hot_encode(sequence)
    return encoded[np.newaxis, ...]
    
def process_fasta_to_batch(fasta_path):
    encoded_list = []
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        
        encoded_seq = prepare_input(seq)
        
        if encoded_seq is not None:
            encoded_list.append(encoded_seq)
    
    if encoded_list:
        batch_array = np.concatenate(encoded_list, axis=0)
        return batch_array
    else:
        return None

enformer_model = hub.load("https://www.kaggle.com/models/deepmind/enformer/TensorFlow2/enformer/1").model
ENFORMER_INPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/enformer_inputs/")
ENFORMER_OUTPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/enformer_outputs/")
target_track = 205

# ['neg, 'log_pos', 'log_test3', 'log_test1', 'log_test2']
for subdir_name in os.listdir(ENFORMER_INPUTS):
    subdir_path = ENFORMER_INPUTS / subdir_name
    if not subdir_path.is_dir():
        continue

    for fasta_path in subdir_path.glob("*.fa"):
        print(f"--- Processing: {fasta_path} ---", flush=True)

        # output paths
        stem = fasta_path.stem.replace("_500bp", "")
        out_subdir = ENFORMER_OUTPUTS / subdir_name
        out_subdir.mkdir(parents=True, exist_ok=True)
        
        save_path = out_subdir / f"{stem}_enformer_{target_track}_preds.npy"

        # check if output already exists
        if save_path.exists() and save_path.stat().st_size > 0:
            print(f"Skipping: {save_path.name} already exists.", flush=True)
            continue

        results = [] 
        x = process_fasta_to_batch(str(fasta_path)) 
        
        for i in range(len(x)):
            # eval one sequence at a time: Shape (1, 196608, 4)
            single_input = x[i : i+1]
            
            pred = enformer_model.predict_on_batch(single_input)
            
            track_subset = pred['mouse'][0, :, target_track]
            
            results.append(track_subset.numpy() if hasattr(track_subset, 'numpy') else track_subset)
            
            if i % 10 == 0:
                print(f"  Sequence {i}/{len(x)}", flush=True)

        # shape (Num_Peaks, 896)
        final_output = np.array(results)
        
        np.save(save_path, final_output)

