# Enformer Analysis

This document outlines how to run the Enformer analysis pipeline.

---

## Overview

The pipeline consists of:

1. Creating the conda environment 
2. Running Enformer predictions
3. Computing correlations on model outputs

---

## File Descriptions

- `enformer_config.yaml`  
  Conda env configuration file

- `run_enformer.py`  
  Main script for running Enformer predictions.

- `submit_enformer.sh`  
  SLURM submission script that wraps `run_enformer.py`.

- `correlations.py`  
  Computes correlations between Enformer predictions and target data.

---

## 1. Create the enformer conda environment

Create the conda env from the config file:

```bash
conda env create -f enformer_config.yaml
conda activate enformer_env
```

## 2. Run Enformer Predictions

Submit the SLURM job:

```bash
sbatch submit_enformer.sh
```

This scripts calls:

```bash
run_enformer.py
```

Make sure to update the paths in `run_enformer.py`

## 3. Compute Correlations

After predictions are complete, run:

```bash
python correlations.py
```


