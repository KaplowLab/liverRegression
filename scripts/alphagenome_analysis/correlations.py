"""
Correlation and Fold-Change Analysis on AlphaGenome predictions vs. real signal values
---------------------------------------------------
This script evaluates AlphaGenome model performance by comparing predicted
chromatin accessibility against ground-truth experimental data across multiple 
species (Mouse, Macaque, Rat, Cow, Pig).

Core Functionalities:
1.  Global Correlation: Calculates Pearson and Spearman coefficients between 
    experimental BED signal and model-predicted .npy values.
2.  Descriptive Statistics: Computes global means, variances, and counts for
    predictions across different sub-cohorts.
3.  Fold-Change (Species-to-Species): Analyzes the model's ability to capture
    evolutionary differences by correlating predicted vs. true log-fold changes 
    between Mouse and other species using one-to-one orthologous peaks.

Outputs:
- correlations.tsv: Statistical metrics per subdir and species.
- all_metadata.tsv: Distribution stats (mean/std) for the predictions.
- species_foldchange.tsv: Mouse-species foldchange analysis.
"""

import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import pearsonr, spearmanr
from pathlib import Path
from sklearn.metrics import mean_squared_error
import os

# PATH DECLARATIONS
INPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/alphagenome_inputs/")
OUTPUTS = Path("/home/azstephe/liverRegression/regression_liver/data/alphagenome_outputs/")
META_DIR = Path("/home/azstephe/liverRegression/regression_liver/data/alphagenome_meta/")
ONE_TO_ONE_DIR = Path("/home/azstephe/liverRegression/regression_liver/data/test_splits/oneToOnePeaks/")

META_DIR.mkdir(parents=True, exist_ok=True)

results_list = []
meta_results = []
species_list = ['macaque', 'rat', 'cow', 'pig']

# HELPER FUNCTIONS
def format_result_value(metric, value):
    if "Count" in metric:
        return int(value)
    elif "P-Val" in metric:
        return f"{value:.2E}"
    elif "%" in metric or metric in ["Pearson", "Spearman", "MSE"]:
        return round(float(value), 3)
    return value

def load_processed_npy(path):
    """Loads the npy and ensures it is a flat 1D array."""
    data = np.load(path)
    return data.flatten()

# GLOBAL CORRELATION & METADATA
# Extracts mean & variance of AlphaGenome predictions
# Extract pearson and spearman correlation coefficients

for subdir_name in os.listdir(OUTPUTS):
    subdir_path = OUTPUTS / subdir_name
    if not subdir_path.is_dir():
        continue

    meta_collector = {}

    for pred_path in subdir_path.glob("*.npy"):
        species = pred_path.name.split('_')[0]
        
        # Load the pre-averaged data
        preds = load_processed_npy(pred_path)

        if species not in meta_collector:
            meta_collector[species] = []
        meta_collector[species].extend(preds)

        # Run on positive groups only
	if subdir_name not in ['log_test1', 'neg']:
            bed_path = INPUTS / subdir_name / f'{species}_liver_TEST_1Mb.bed'

            if bed_path.exists():
                true_values = pd.read_csv(bed_path, header=None, sep=r'\s+').iloc[:, 4].values

                if len(preds) != len(true_values):
                    print(f"Warning: Size mismatch in {subdir_name}! Preds: {len(preds)}, True: {len(true_values)}")
                    min_len = min(len(preds), len(true_values))
                    preds = preds[:min_len]
                    true_values = true_values[:min_len]

                r, p = pearsonr(true_values, preds)
                s, sp = spearmanr(true_values, preds)
                results_list.append({
                    'subdir': subdir_name, 'species': species,
                    'pearson_r': r, 'pearson_p': p,
                    'spearman_rho': s, 'spearman_p': sp,
                    'n': len(preds)
                })

    for species, all_preds in meta_collector.items():
        if all_preds:
            meta_results.append({
                'subdir': subdir_name,
                'species': species,
                'global_avg': np.mean(all_preds),
                'std_dev': np.std(all_preds),
                'count': len(all_preds)
            })

pd.DataFrame(results_list).to_csv(META_DIR / "correlations.tsv", sep='\t', index=False)
pd.DataFrame(meta_results).to_csv(META_DIR / "all_metadata.tsv", sep='\t', index=False)


# FOLDCHANGE ANALYSIS

def get_signal_dict(bed_path, pred_path):
    bed = pd.read_csv(bed_path, sep=r'\s+', header=None)
    true_map = dict(zip(bed[3], bed[4]))
    
    # Using the pre-averaged npy
    preds = load_processed_npy(pred_path)
    pred_map = dict(zip(bed[3], preds))
    
    return true_map, pred_map

rows = []
for species in species_list:
    bridge_path = ONE_TO_ONE_DIR / f"{species}_mouse.bed"
    if not bridge_path.exists(): continue
    bridge = pd.read_csv(bridge_path, header=None, sep='\t')

    m_true, m_pred = get_signal_dict(
        INPUTS / "log_pos" / "mouse_liver_TEST_1Mb.bed",
        OUTPUTS / "log_pos" / "mouse_liver_TEST_human_0002107_preds.npy" 
    )

    s_true, s_pred = get_signal_dict(
        INPUTS / "log_test2" / f"{species}_liver_TEST_1Mb.bed",
        OUTPUTS / "log_test2" / f"{species}_liver_TEST_human_0002107_preds.npy"
    )

    data = []
    for _, row in bridge.iterrows():
        s_id, m_id = row[4], row[14]
        if m_id in m_true and s_id in s_true:
            data.append({
                'm_true': m_true[m_id], 'm_pred': m_pred[m_id],
                's_true': s_true[s_id], 's_pred': s_pred[s_id]
            })

    df = pd.DataFrame(data)
    if df.empty: continue

    # Log-fold change calculation 
    x = df['m_true'] - df['s_true']
    y = df['m_pred'] - df['s_pred']

    pearson, pp = pearsonr(x, y)
    spearman, ps = spearmanr(x, y)
    mse = mean_squared_error(x, y)

    num_ss = (np.sign(x) == np.sign(y)).sum()
    perc_ss = num_ss / len(df)

    metrics = [
        ('Same Sign Count', num_ss),
        ('Total Count', len(df)),
        ('Same Sign %', perc_ss),
        ('Pearson', pearson),
        ('Pearson P-Val', pp),
        ('Spearman', spearman),
        ('Spearman P-Val', ps),
        ('MSE', mse)
    ]

    for metric_name, value in metrics:
        rows.append({
            'Species': species,
            'Group': 'Test2',
            'Metric': metric_name,
            'AlphaGenome': format_result_value(metric_name, value)
        })

pd.DataFrame(rows).to_csv(META_DIR / "foldchange_results.tsv", sep='\t', index=False)
