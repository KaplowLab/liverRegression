# Motif Discovery (MEME-ChIP)
[MEME Suite](https://meme-suite.org/meme/) was used to do motif discovery.

This section provides example code as a template.

### Analysis Examples:

**Enrichment**: Running MEME-ChIP on the strongest mouse OCRs.

**Differential**: Comparing the strongest peaks against a background of weak peaks.

#### Filter for only strong peaks with `OCROrthologPrediction`
``` bash
# Assumes a CSV file of the strongest peak names is already generated.
python ~/repos/OCROrthologPrediction/src/filterPeakName.py \
    --unfilteredPeakFileName ~/data/raw/mouse_liver_pos_ALL.narrowPeak \
    --peakListFileName ~/data/mouse_certainty/mouse_strong_peaklist.csv \
    --peakNameCol 0 \
    --outputFileName ~/data/mouse_certainty/mouse_liver_pos_STRONG.narrowPeak
```
#### Summit center and expand to 500 bp with `atac_data_pipeline`
``` bash
python ~/repos/atac_data_pipeline/scripts/preprocessing.py expand_peaks \
-l 500 \
-i ~/data/mouse_certainty/mouse_liver_pos_STRONG.narrowPeak \
-o ~/data/mouse_certainty/mouse_liver_pos_STRONG_500bp.narrowPeak
```
#### Extract sequences with `bedtools`
``` bash
bedtools getfasta -fi /~/data/MouseGenome/mm10.fa \
    -bed ~/data/mouse_certainty/mouse_liver_pos_STRONG_500bp.narrowPeak \
    > ~/data/mouse_certainty/mouse_liver_pos_STRONG_500bp.fa
```

#### Run MEME-ChIP
Custom sbatch script provided in `scripts/memechip.sh`

Includes both direct shell commands and SLURM (sbatch) script examples for cluster submission.
``` bash
# Sbatch command for memechip on strong
sbatch --mem=4G -o strong.o -J strong -p gpu -n 1 --gres gpu:1 ~/scripts/memechip.sh \
    ~/data/memechip/mouseStrongMemechip \
    ~/data/mouse_certainty/mouse_liver_pos_STRONG_500bp.fa

# Strong vs. weak meme-chip run
~/anaconda2/bin/meme-chip -o ~/data/mouse_strength/mouseStrongVsWeakMemechip \
    -db ~/data/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme \
    -neg ~/data/mouse_certainty/mouse_liver_pos_WEAK_500bp.fa \
    ~/data/mouse_certainty/mouse_liver_pos_STRONG_500bp.fa
```
