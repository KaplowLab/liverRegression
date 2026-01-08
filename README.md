# LiverRegression
This repository contains the codebase and data processing pipelines used for the project: "Challenges in Predicting Chromatin Accessibility Differences between Species."

The project focuses on using regression models to predict species-specific open chromatin regions, utilizing specialized normalization techniques to handle cross-species signal variation.

## Key Features
Quantile Normalization (quantile_normalize.py): Standardizes signal distributions across species for fixed-size peak sets.

Extended Quantile Normalization (EQN.py): A robust normalization strategy that handles varying peak counts across different species.

## Data Preprocessing
Detailed documentation on data generation can be found in the `workflows/` directory. This includes the step-by-step logic for log, QN, EQN-transformation, summit-centering, and base-pair expansion (500bp/2000bp).

### Data Split Organization
Processed data is organized into Train, Validation, and Test splits. Below is an example summary of the directory structure within ~/data/:

| Split Type | Subset | Description |
| :--- | :--- | :--- |
| **Train** | `train_splits/` | Includes `neg` (negative sets) and `log_pos` (log-transformed positives). |
| **Validation** | `val_splits/` | Includes `neg`, `log_pos`, and subsets (`val1`, `log_val2`, `log_val3`). |
| **Test** | `test_splits/` | Includes `neg`, `log_pos`, and subsets (`test1`, `log_test2`, `log_test3`). |

## Training Regression Models
[Pfenning Lab CNN Pipeline](https://github.com/pfenninglab/cnn_pipeline) was used for model training. 

Example configuration and hyperparameter sweep files are located in the `examples/` directory.

#### Example Usage
```
# Single training
bash train.sh examples/config.yaml

# Hyperparameter sweep
bash start_sweep.sh examples/sweep-config.yaml
bash start_agents.sh 100 4 <sweep_id>
```
Models and logs are saved to: `~/repos/cnn_pipeline/wandb/`

An example run output is provided in: `examples/run-20250225_222915-bdbi7l3n`

## Model Activations
`scripts/get_activations.sh` automates the extraction of model predictions for all datasets. This script automatically submits sbatch jobs to the cluster and uses a GPU. Gets activations for the following set: 

Positive and Negative sets for Training, Validation, and Testing.
Subsets: Val 1, 2, 3 and Test 1, 2, 3.

Note on Usage: These scripts use project-specific naming conventions and data structures. While they serve as documented templates for the analysis methodology, they will require adaptation to work with external datasets.
```
bash get_activations.sh run-20250225_222915-bdbi7l3n
```
Output Directory: ~/data/model_outputs/${ID}_FINAL

## Model Evaluations

Examples for correlating model predictions with observed signal values are available in the provided Python notebooks.

Note on Usage: These scripts use project-specific naming conventions and data structures. While they serve as documented templates for the analysis methodology, they will require adaptation to work with external datasets.

`scripts/foldchange_eval.ipynb`
`scripts/prediction_eval.ipynb`

## Memechip Analysis
