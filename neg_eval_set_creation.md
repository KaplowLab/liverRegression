# Workflow: Generation of Negative Evaluation Sets for Cross-Species OCR Prediction
This documentation outlines the bioinformatic pipeline used to create the negative (non-enhancer) evaluation sets. This process identifies the rat regions that are inaccessible, whose mouse orthologs are accessible.
The same pipeline is applied to generate training and validation sets; simply modify the chromosomal partitions in Step 3 accordingly.

## 1. Cross-Species Peak Mapping (Mouse to Rat Liftover)
We use HALPER to map Mouse peak coordinates to the Rat genome. This identifies the orthologous genomic locations in Rat for all known Mouse OCRs.
```bash
# Map all mouse peaks to rat using HALPER via SLURM
sbatch -p pfen1 -w compute-1-11 -o ~/mouse_rat.o \
  /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -s Mus_musculus \
  -t Rattus_norvegicus \
  -o /scratch/azstephe/mouseToRat/ \
  -b /home/azstephe/liverRegression/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak \
  --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover
```
Primary output: mouse_liver_pos_ALL.Mus_musculusToRattus_norvegicus.HALPER.narrowPeak.gz

## 2. Create Test Set 1
We identify regions that are accessible in Mouse but are not accessible in Rat.

### Identify inaccessible Rat region coordinates whose orthologs are accessible in mouse
We achieve this by taking the Mouse peaks mapped to the Rat genome and subtracting all Rat OCRs.

```bash
# Subtract known rat peaks from mapped mouse peaks
# -A ensures any overlap results in the removal of the entire peak
bedtools subtract -A \
  -a /home/azstephe/liverRegression/regression_liver/data/mapped/mouse_liver_pos_ALL.Mus_musculusToRattus_norvegicus.HALPER.narrowPeak.gz \
  -b /home/azstephe/regression_liver/data/raw/rat_liver_pos_ALL.narrowPeak \
  > /home/azstephe/liverRegression/regression_liver/data/mapped/mouseToRat_liver_ratNon_mouseEnhancer_ALL.narrowPeak
```
Output: A list of Mouse peaks that are non-enhancers in Rat.

### Filter for Test Set
Filter the previous output to keep only the coordinates whose Mouse peaks belong to the test set (Mouse chr1 and chr2).
```bash 
# Filter for Rat coordinates whose mouse orthologs are in the Test Set
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py \
  --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/mapped/mouseToRat_liver_ratNon_mouseEnhancer_ALL.narrowPeak \
  --peakListFileName /home/azstephe/liverRegression/regression_liver/data/raw/mouse_liver_pos_TEST.narrowPeak \
  --peakNameCol 3 \
  --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test1/rat_liver_TEST.narrowPeak
```
### Peak Standardization
Center the identified peaks on their summits and expand them to 500bp to prepare them for the regression model. Performed using the Pfenning Lab `atac_data_pipeline` repository.

```bash
# Center and expand peaks to 500bp for model input
python /home/azstephe/repos/atac_data_pipeline/scripts/preprocessing.py expand_peaks \
  -l 500 \
  -i /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test1/rat_liver_TEST.narrowPeak \
  -o /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test1/rat_liver_TEST_500bp.bed
```
Primary output used in training: /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test1/rat_liver_TEST_500bp.bed
