import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

ENFORMER_OUT = Path('/home/azstephe/liverRegression/regression_liver/data/enformer_outputs')
ENFORMER_INPUTS = Path('/home/azstephe/liverRegression/regression_liver/data/enformer_inputs/')
target_track = 205

results_list = []
neg_results = []

dir_list = ['log_test2', 'log_test1']

for subdir_name in dir_list:
    subdir_path = ENFORMER_OUT / subdir_name
    if not subdir_path.is_dir(): 
        continue

    # GET METADATA FOR ALL SUBDIRS
    species_neg_collector = {}
        
    for pred_path in subdir_path.glob("*_205_preds.npy"):
        # Extract species (e.g., 'human' or 'mouse') from the filename
        species = pred_path.name.split('_')[0]
            
        preds_raw = np.load(pred_path) # shape (N, 896)
        # Average the center 4 bins for each sequence in this file
        avg_per_seq = np.mean(preds_raw[:, 446:450], axis=1)
            
        if species not in species_neg_collector:
            species_neg_collector[species] = []
            
        species_neg_collector[species].extend(avg_per_seq)
        
    for species, all_preds in species_neg_collector.items():
        if all_preds:
            neg_results.append({
                'subdir': subdir_name,
                'species': species,
                'global_avg': np.mean(all_preds),
                'std_dev': np.std(all_preds),
                'count': len(all_preds)
            })
            print(f"Processed {subdir_name} - {species} as Negative Control.", flush=True)

    # DO THIS FOR THE POSITIVES
    if subdir_name not in ['log_test1', 'neg']:
        for pred_path in subdir_path.glob("*_205_preds.npy"):
            species = pred_path.name.split('_')[0]
            preds_raw = np.load(pred_path)
            preds_avg = np.mean(preds_raw[:, 446:450], axis=1).flatten()

            bed_path = ENFORMER_INPUTS / subdir_name / f'{species}_liver_TEST_500bp.bed'
            if not bed_path.exists(): continue
            
            true_values = pd.read_csv(bed_path, header=None, sep=r'\s+').iloc[:, 4].values

            if len(preds_avg) != len(true_values):
                min_len = min(len(preds_avg), len(true_values))
                preds_avg, true_values = preds_avg[:min_len], true_values[:min_len]
            
            if len(preds_avg) > 1:
                r, p = pearsonr(true_values, preds_avg)
                s, sp = spearmanr(true_values, preds_avg)
                results_list.append({
                    'subdir': subdir_name, 'species': species,
                    'pearson_r': r, 'pearson_p': p, 'spearman_rho': s, 'spearman_p': sp, 'n': len(preds_avg)
                })

pd.DataFrame(results_list).to_csv("/home/azstephe/liverRegression/regression_liver/data/enformer_meta/correlations.tsv", sep='\t', index=False)
pd.DataFrame(neg_results).to_csv("/home/azstephe/liverRegression/regression_liver/data/enformer_meta/all_metadata.tsv", sep='\t', index=False)
