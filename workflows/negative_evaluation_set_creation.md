# Workflow: Generation of Negative Evaluation Sets for Cross-Species OCR Prediction

This documentation outlines the bioinformatics pipeline used to generate the negative evaluation sets from ATAC-seq data. 

In this example, mouse (*Mus musculus*) serves as the training species, while rat (*Rattus norvegicus*) serves as the evaluation species.

Below is the step-by-step implementation for generating the negative evaluation partitions: Neg Test and Test 1.

The same pipeline is applied to generate training and validation sets; simply modify the chromosomal partitions in the filtering of Step 2 accordingly.

## Create Negative Test Set

The creation of the negative test set is described in [Kaplow et al.](https://www.science.org/doi/10.1126/science.abm7993), *Science* (2023). Here we will briefly describe the process. 
### 1. Map candidate enhancers from all non-rat species to rat
Using the `halLiftover-postprocessing` repository, map mouse peak coordinates to the rat genome. Do this for all non-rat species (in our case: mouse, macaque, cow, and pig)
```bash
# Map all mouse peaks to rat using HALPER via SLURM
sbatch -p [partition] -w [node] -o ~PROJECT_DIR/mouseToRat.o \
  $PROJECT_DIR/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -s Mus_musculus \
  -t Rattus_norvegicus \
  -o $PROJECT_DIR/data/mouseToRat/ \
  -b $PROJECT_DIR/data/candidate_enhancers/mouse_liver_pos_ALL.narrowPeak
```
Primary output: $PROJECT_DIR/data/mouseToRat/mouse_liver_pos_ALL.Mus_musculusToRattus_norvegicus.HALPER.narrowPeak.gz

### 2. Remove those that overlap with rat OCRs (both reproducible and nonreproducible)
### 3. Union the resulting regions
### 4. Split according to the train-val-test split

Primary output used in training: $PROJECT_DIR/data/test_splits/neg/rat_liver_TEST_500bp.bed

## Create Test Set 1
Regions that are accessible in mouse but are not accessible in rat.

### 1. Identify inaccessible rat region coordinates whose orthologs are accessible in mouse
Take the mouse peaks mapped to the rat genome and subtracting all rat OCRs (reproducible and non-reproducible).
```bash
# Subtract known rat peaks from mapped mouse peaks
# -A ensures any overlap results in the removal of the entire peak
bedtools subtract -A \
  -a $PROJECT_DIR/data/mouseToRat/mouse_liver_pos_ALL.Mus_musculusToRattus_norvegicus.HALPER.narrowPeak.gz \
  -b $PROJECT_DIR/data/candidate_enhancers/rat_liver_pos_ALL.narrowPeak \
  > $PROJECT_DIR/data/mouseToRat/mouseToRat_liver_ratNon_mouseEnhancer_ALL.narrowPeak
```
Output: A list of mouse peaks that are non-enhancers in rat.

### 2. Filter for Test Set
Filter the previous output to keep only the coordinates whose mouse peaks belong to the test set (mouse chr1 and chr2).
```bash 
# Filter for rat coordinates whose mouse orthologs are in the Test Set
python $PROJECT_DIR/repos/OCROrthologPrediction/src/filterPeakName.py \
  --unfilteredPeakFileName $PROJECT_DIR/data/mouseToRat/mouseToRat_liver_ratNon_mouseEnhancer_ALL.narrowPeak \
  --peakListFileName $PROJECT_DIR/data/candidate_enhancers/mouse_liver_pos_TEST.narrowPeak \
  --peakNameCol 3 \
  --outputFileName $PROJECT_DIR/data/test_splits/log_test1/rat_liver_TEST.narrowPeak
```
### Peak Standardization
Using the `atac_data_pipeline` repository, center the identified peaks on their summits and expand them to 500bp to prepare them for the regression model. 

```bash
# Center and expand peaks to 500bp for model input
python $PROJECT_DIR/repos/atac_data_pipeline/scripts/preprocessing.py expand_peaks \
  -l 500 \
  -i $PROJECT_DIR/data/test_splits/log_test1/rat_liver_TEST.narrowPeak \
  -o $PROJECT_DIR/data/test_splits/log_test1/rat_liver_TEST_500bp.bed
```
Primary output used in training: $PROJECT_DIR/data/test_splits/log_test1/rat_liver_TEST_500bp.bed
