"""
Predict with AlphaGenome
------------------------------------------
This script passes genomic sequences (FASTA) into the AlphaGenome DNA model
to predict chromatin accessibility (ATAC) specifically for liver tissue 
(UBERON:0002107) in human.

Workflow:
1. Iterates through predefined subdirectories in the input path.
2. Loads .fa files and pads sequences to the required 1Mb length.
3. Queries the AlphaGenome API for ATAC-seq predictions in Homo Sapiens.
4. Extracts the central 500bp signal from the model output.
5. Calculates the mean prediction and saves the results as NumPy arrays (.npy).
"""

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
import os
import numpy as np
from Bio import SeqIO
from pathlib import Path
import time
from datetime import timedelta

# CONFIGURATION
ALPHA_INPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/alphagenome_inputs")
ALPHA_OUTPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/alphagenome_outputs/")

API_KEY = 'AIzaSyBHUe7xJ3qMQAM-zFuxC5ESQttfQ6TaSLI'

dna_model = dna_client.create(API_KEY)

for subdir_name in ['log_test3', 'log_pos', 'neg', 'log_test1', 'log_test2']:

    subdir_path = ALPHA_INPUTS / subdir_name
    if not subdir_path.is_dir(): 
        continue

    for fasta_path in subdir_path.glob("*.fa"):
        stem = fasta_path.stem.replace("_1Mb", "")
        out_subdir = ALPHA_OUTPUTS / subdir_name
        save_path = out_subdir / f"{stem}_human_0002107_preds.npy"
        
        if save_path.exists():
            print(f"Skipping: {save_path} already exists.")
            continue
            
        print(f"\n Processing: {fasta_path.name}")
        results = []
        
        total_seqs = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))

        for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
            raw_seq = str(record.seq).upper()
            padded_seq = raw_seq.center(dna_client.SEQUENCE_LENGTH_1MB, 'N')
            
            try:
                output = dna_model.predict_sequence(
                    sequence=padded_seq,
                    requested_outputs=[dna_client.OutputType.ATAC],
                    organism=dna_client.Organism.HOMO_SAPIENS,
                    ontology_terms=['UBERON:0002107'], 
                )

                preds = output.atac.values 
                total_len = preds.shape[0]
                
                # Get center 500bp
                center = total_len // 2
                start = center - 250
                end = center + 250
                
                center_signal = preds[start:end, :]
                avg_pred = center_signal.mean()
                results.append(avg_pred)

            except Exception as e:
                print(f"   Error processing sequence {i}: {e}")
                continue
            
            # Progress and ETA reporting every 50 sequences
            if i > 0 and i % 50 == 0:
                print(f"Progress: {i}/{total_seqs}", flush=True)

        # Save results
        if results:
            out_subdir.mkdir(parents=True, exist_ok=True)
            np.save(save_path, np.array(results))
            
            print(f"Saved {len(results)} predictions to: {save_path}")
        else:
            print(f"No results generated for {fasta_path}")

print("\n All AlphaGenome inferences complete.")
