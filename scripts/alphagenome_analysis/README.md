# AlphaGenome Analysis

This document outlines the steps required to run the AlphaGenome analysis pipeline.

## Installation Notes

Installing AlphaGenome required extensive troubleshooting. The setup described below works, but JAX failed to load with GPU support, so all computations were performed on CPUs instead of GPUs.

---

## 1. Build the Conda Environment

Create and activate the environment:

```bash
bash build_alpha_env.sh
conda activate alphagenome_env
```

## 2. Run AlphaGenome Predictions

Submit the SLURM job:
```bash
sbatch submit_alphagenome.sh
```

This script is a wrapper for:
```bash
run_alphagenome.py
```

## 3. Compute correlations

After predictions are complete, calculate correlations:
```bash
python correlations.py
```

