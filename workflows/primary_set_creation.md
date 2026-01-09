# Workflow for Primary File Creation
This document outlines the bioinformatics pipeline used to generate and preprocess primary data files. The methods assume the files in `processed/` are the processed datasets as described in [Kaplow et al.](https://link.springer.com/article/10.1186/s12864-022-08450-7) (2022) ("Inferring mammalian tissue-specific regulatory conservation by predicting tissue-specific differences in open chromatin").

All operations are assumed to be executed within the data/ directory.

## NarrowPeak File Format
| Column | Name | Type | Description |
| :--- | :--- | :--- | :--- |
| **1** | `chrom` | string | Name of the chromosome (e.g., chr1, chrX). |
| **2** | `chromStart` | int | The starting position of the peak (0-based). |
| **3** | `chromEnd` | int | The ending position of the peak (1-based). |
| **4** | `name` | string | Name of the peak (e.g., `peak1`). |
| **5** | `score` | int | Indicates how dark the peak will be displayed in a browser (0–1000). |
| **6** | `strand` | char | Orientation (+ or -). Use `.` if no strand is assigned. |
| **7** | `signalValue` | float | Measurement of enrichment (e.g., fold-change). **Used for EQN.** |
| **8** | `pValue` | float | Statistical significance as $-log_{10}(p-value)$. |
| **9** | `qValue` | float | Statistical significance (FDR) as $-log_{10}(q-value)$. |
| **10** | `peak` | int | **Summit offset:** 0-based offset from `chromStart` to the peak summit. |

We are most interested in Columns 1, 2, 3, 4, 7.

## Signal Preprocessing Strategies
Processed `narrowPeak` data must be transformed to stabilize variance and ensure comparability across species. Multiple strategies are explored in this paper and helper functions are provided in `scripts/`. Here, we give examples of the transformations to normalized signal values.

### 1. Log Transformation
This step applies a $log(x+1)$ transformation to the signal values (Column 7 in narrowPeak format). This preserves the original file structure.
``` bash
# Example: Logging signal values
python scripts/log.py ~/data/processed/mouse_liver_pos_ALL.narrowPeak
```
Output Directory: ~/data/log

### 2. Quantile Normalization (QN)
To perform standard Quantile Normalization, all input files must have the same number of peaks. Various straetegies to downsample larger datasets to match the species with the smallest peak count include random, distribution matched random, and keeping the strongest signals. 

Note: In this specific dataset, the Pig sample contained the fewest peaks (20,615). 20,615 peaks were randomly selected from all other species and stored in data/log_same_size/.

``` bash
# Example: Running QN across five species. Must include all the species files at once.
python scripts/quantile_normalize.py \
  ~/data/log/mouse_liver_pos_ALL.narrowPeak \
  ~/data/log/macaque_liver_pos_ALL.narrowPeak \
  ~/data/log/rat_liver_pos_ALL.narrowPeak \
  ~/data/log/cow_liver_pos_ALL.narrowPeak \
  ~/data/log/pig_liver_pos_ALL.narrowPeak
```
Output Directory: `~/data/quantile_norm/`

### 3. Extended Quantile Normalization (EQN)
Unlike standard QN, Extended Quantile Normalization allows for peak files of varying lengths.

Prerequisites: This script assumes the existence of log/, log_sorted/, and quantile_norm/ subdirectories.

``` bash
# Example: Applying EQN to the Mouse dataset. Do this for all species.
python scripts/EQN.py mouse_liver_pos_ALL.bed ~/data/
```

Output Directory: `~/data/eqn/`

## Genomic Coordinate Processing
Once signals are normalized, the OCR lengths must be standardized for model input.

### Base Pair Extension and Formatting
Center each peak on its summit and expand it to a uniform width using the Pfenning Lab atac_data_pipeline. The lengths 500bp and 2000bp were used in the paper. Do this for all species.

``` bash
# 1. Summit center and expand to 500bp
python ~/repos/atac_data_pipeline/scripts/preprocessing.py expand_peaks \
  -l 500 \
  -i ~/data/log/mouse_liver_pos_ALL.narrowPeak \
  -o ~/data/log/mouse_liver_pos_ALL_500bp.narrowPeak 

# 2. Extract essential columns (BED5 format: chr, start, stop, name, signal)
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $7}' \
  ~/data/log/mouse_liver_pos_ALL_500bp.narrowPeak > \
  ~/data/log/mouse_liver_pos_ALL_500bp.bed
```
Primary output used for evaluation set creation: `~/data/log/mouse_liver_pos_ALL_500bp.bed`
