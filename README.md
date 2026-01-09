# LiverRegression
This repository contains the codebase and data processing pipelines used for the project: 

### Challenges in Predicting Chromatin Accessibility Differences between Species

Amy Stephen, Arian Raje, Heather H. Sestili, Morgan E. Wirthlin, Alyssa J. Lawler, Ashley R. Brown, William R. Stauffer, Andreas R. Pfenning, Irene M. Kaplow. bioRxiv 2025; doi: https://doi.org/10.1101/2025.11.09.687449

https://www.biorxiv.org/content/10.1101/2025.11.09.687449v1

The project highlights the challenges of using regression convolutional neural networks (CNNs) to predict chromatin accessibility differences between species in bulk liver ATAC-seq data.

Note on Usage: Most scripts use project-specific naming conventions and data structures. While they serve as documented templates for the analysis methodology, they will require adaptation to work with external datasets.

## Key Features
Quantile Normalization (quantile_normalize.py): Standardizes signal distributions across species for fixed-size peak sets.

Extended Quantile Normalization (EQN.py): A robust normalization strategy that handles varying peak counts across different species.

## External Tools

A number of public tools were used in this work. They are linked below.

[atac_data_pipeline](https://github.com/pfenninglab/atac_data_pipeline)

[halLiftover-postprocessing](https://github.com/pfenninglab/halLiftover-postprocessing)

[OCROrthologPrediction](https://github.com/pfenninglab/OCROrthologPrediction)

[cnn_pipeline](https://github.com/pfenninglab/cnn_pipeline)

[MEME-Suite](https://meme-suite.org/meme/tools)

## Data Preprocessing
Detailed documentation on data generation can be found in the `workflows/` directory. This includes the step-by-step logic for log, QN, EQN-transformation, summit-centering, and base-pair expansion (500bp/2000bp).

### Data Split Organization
Processed data is organized into Train, Validation, and Test splits. Below is an example summary of the directory structure within `$PROJECT_DIR/data/`:

| Split Type | Subset | Description |
| :--- | :--- | :--- |
| **Train** | `train_splits/` | Includes `neg` (negative sets) and `log_pos` (log-transformed positives). |
| **Validation** | `val_splits/` | Includes `neg`, `log_pos`, and subsets (`val1`, `log_val2`, `log_val3`). |
| **Test** | `test_splits/` | Includes `neg`, `log_pos`, and subsets (`test1`, `log_test2`, `log_test3`). |

## Training Regression Models
The [Pfenning Lab CNN Pipeline](https://github.com/pfenninglab/cnn_pipeline) was used for model training. 

Example configuration and hyperparameter sweep files are located in the `examples/` directory.

#### Example Usage
```
# Single training
bash train.sh examples/config.yaml

# Hyperparameter sweep
bash start_sweep.sh examples/sweep-config.yaml
bash start_agents.sh 100 4 <sweep_id>
```
Models and logs are saved to: `$PROJECT_DIR/repos/cnn_pipeline/wandb/`

An example run output follows the following file structure:
```
run-20250225_222915-bdbi7l3n
├── files
│   ├── conda-environment.yaml
│   ├── config.yaml
│   ├── media
│   │   └── graph
│   │       └── graph_summary_c62f178e911f2edbadcc.graph.json
│   ├── model-best.h5
│   ├── model-latest.h5
│   ├── output.log
│   ├── requirements.txt
│   ├── wandb-metadata.json
│   └── wandb-summary.json
├── logs
│   ├── debug-internal.log
│   └── debug.log
├── run-bdbi7l3n.wandb
└── tmp
    └── code
```

## Model Activations
`scripts/get_activations.sh` automates the extraction of model predictions for all datasets. This script automatically submits sbatch jobs to the cluster and uses a GPU. Gets activations for the following set: 

Positive and Negative sets for Training, Validation, and Testing.
Subsets: Val 1, 2, 3 and Test 1, 2, 3.

```
bash get_activations.sh run-20250225_222915-bdbi7l3n
```
Output Directory: `$PROJECT_DIR/data/model_outputs/${ID}_FINAL` where `$ID=bdbi7l3n`

## Model Evaluations

Examples for correlating model predictions with observed signal values are available in the provided Python notebooks.

`notebooks/foldchange_eval.ipynb`
`notebooks/prediction_eval.ipynb`

## Motif Enrichment Analysis
Motif Enrichment was conducted using [MEME-ChIP](https://meme-suite.org/meme/tools/meme-chip). 

Example code is provided in `workflow/mouse_memechip.md`.

The Motif Enrichment results for the project can be found [here](http://daphne.compbio.cs.cmu.edu/files/azstephe/liver_regression_resource/)

## Dependencies
Dependencies are consistent with those of the Pfenning Lab [CNN Pipeline](https://github.com/pfenninglab/cnn_pipeline). 

## Contact

Irene Kaplow (ikaplow@cs.cmu.edu)

Andreas Pfenning (apfenning@cmu.edu)

Amy Stephen (astephen@berkeley.edu)
