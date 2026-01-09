# Workflow: Generation of Positive Evaluation Sets for Cross-Species OCR Prediction

This documentation outlines the bioinformatics pipeline used to generate positive evaluation sets from ATAC-seq data. 

In this example, Mouse (*Mus musculus*) serves as the training species, while Rat (*Rattus norvegicus*) serves as the evaluation species.

Below is the step-by-step implementation for generating the positive test evaluation partitions: Test, Test 2, and Test 3.

The same pipeline is applied to generate training and validation sets; simply modify the chromosomal partitions in Step 2 accordingly.

## 1. Cross-Species Peak Mapping (Liftover)
Using `HALPER`, map rat peak coordinates to the mouse genome. It is important to use the `$PROJECT_DIR/enhancer_candidates/` files before any summit-centering and peak expansion to correctly map enhancer candidate orthologs. The output files are very large and are often stored in a scratch directory which will be omitted for clarity.
Remember to first `conda activate hal`

```bash
# Map all rat peaks to mouse using HALPER via SLURM
sbatch -p pfen1 -w compute-1-11 -o $PROJECT_DIR/ratToMouse.o \
  $PROJECT_DIR/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -s Rattus_norvegicus \
  -t Mus_musculus \
  -o $PROJECT_DIR/data/ratToMouse/ \
  -b $PROJECT_DIR/data/enhancer_candidates/rat_liver_pos_ALL.narrowPeak \
  --halPath ~/src/hal/bin/halLiftover
```

Primary output: `$PROJECT_DIR/data/ratToMouse/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz`

Column 3: Mouse coordinates of the mapped rat peak.
Column 4: Original Rat peak name.

## 2. Create Positive Test Set
Rat regions that are orthologs of mouse chr1 and chr2 serve as our test set.
### Identify Postive Test Peak Set
```bash
# Extract peaks mapped to mouse chr1 and chr2
zcat $PROJECT_DIR/data/ratToMouse/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz \
  | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' \
  > $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_TEST.narrowPeak
```

Note: change the `grep` command to `grep -E '^(chr8[[:space:]]|chr9[[:space:]])'` to generate validation set.

### Filter Positive Test Set
Filter the Test peak set against our log-signal file to get a file with the OCR coordinates and respective logged signal values.
```bash
# Filter test set peaks using the log signal values file
python ~/repos/OCROrthologPrediction/src/filterPeakName.py \
  --unfilteredPeakFileName $PROJECT_DIR/data/log_transformed/rat_liver_pos_ALL_500bp.bed \
  --peakListFileName $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_TEST.narrowPeak \
  --peakNameCol 3 \
  --outputFileName $PROJECT_DIR/data/test_splits/log_pos/rat_liver_TEST_500bp.bed
```

Primary output used in training: `$PROJECT_DIR/data/test_splits/log_pos/rat_liver_TEST_500bp.bed`

## 3. Create Test Set 2
We identify regions that are accessible in both Rat and Mouse. 
### Identify Test 2 Peak Set
Use bedtools intersect with specific flags to preserve the metadata from both species.
```bash
# Intersect mapped rat peaks with existing mouse peaks
# -wa: Write the original entry from file A (Rat mapped to Mouse)
# -wb: Write the original entry from file B (Mouse peak)
bedtools intersect -wa -wb \
  -a $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_TEST.narrowPeak \
  -b $PROJECT_DIR/data/candidate_enhancers/mouse_liver_pos_ALL.narrowPeak \
  > $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb.narrowPeak
```

### Filter Test 2 Set
Filter the Test 2 peak set against our log-signal file to get a file with the OCR coordinates and respective logged signal values.
```bash
# Filter the conserved test set against the log signal values file
python ~/repos/OCROrthologPrediction/src/filterPeakName.py \
  --unfilteredPeakFileName $PROJECT_DIR/data/log_transformed/rat_liver_pos_ALL_500bp.bed \
  --peakListFileName $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb.narrowPeak \
  --peakNameCol 3 \
  --outputFileName $PROJECT_DIR/data/test_splits/log_test2/rat_liver_TEST_500bp.bed
```

Primary output used in training: `$PROJECT_DIR/data/test_splits/log_test2/rat_liver_TEST_500bp.bed`

## 4. Create Test Set 3
We identify regions that are accessible in Rat but not in Mouse.
### Identify Test 3 Peak Set
Subtract all reproducible and non-reproducible mouse peaks from our mapped rat peaks.
```bash
# Subtract mouse peaks from mapped rat peaks
# -A ensures that any overlap results in the removal of the entire rat peak
bedtools subtract -A \
  -a $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_TEST.narrowPeak \
  -b $PROJECT_DIR/data/rep_and_nonRep/mouse_rep_nonRep.narrowPeak.gz \
  > $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.narrowPeak
```

### Filter Test 3 Set
Filter the Test 3 peak set against our log-signal file to get a file with the OCR coordinates and respective logged signal values.
```bash
# Filter test 3 set peaks using the log signal values file
python ~/repos/OCROrthologPrediction/src/filterPeakName.py \
  --unfilteredPeakFileName $PROJECT_DIR/data/log_transformed/rat_liver_pos_ALL_500bp.bed \
  --peakListFileName $PROJECT_DIR/data/ratToMouse/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.narrowPeak \
  --peakNameCol 3 \
  --outputFileName $PROJECT_DIR/data/test_splits/log_test3/rat_liver_TEST_500bp.bed
```

Primary output used in training: `$PROJECT_DIR/data/test_splits/log_test3/rat_liver_TEST_500bp.bed`

