# Workflow: Generation of Evaluation Sets for Cross-Species OCR Prediction

This documentation outlines the bioinformatic pipeline used to create evaluation sets from ATAC-SEQ data. The example workflow maps Open Chromatin Region (OCR) peaks from **Rat** (*Rattus norvegicus*) to **Mouse** (*Mus musculus*) and generates the Test, Test 1, Test 2, and Test 3 sets.

The same pipeline is applied to generate training and validation sets; simply modify the chromosomal partitions in Step 2 accordingly.

## 1. Cross-Species Peak Mapping (Liftover)
Using `HALPER`, we map rat peak coordinates to the mouse genome. This allows us to identify orthologous regions across the two species.



```bash
# Map all rat peaks to mouse using HALPER via SLURM
sbatch -p pfen1 -w compute-1-11 -o ~/rat_mouse.o \
  /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -s Rattus_norvegicus \
  -t Mus_musculus \
  -o /scratch/azstephe/liver/ratToMouse/ \
  -b /home/azstephe/liverRegression/regression_liver/data/raw/rat_liver_pos_ALL.narrowPeak \
  --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover
```

Primary output: `rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz`

Column 3: Mouse coordinates of the mapped rat peak.
Column 4: Original Rat peak name.

## 2. Identify Postive Test Peak Set (Mouse Chromosomes 1 & 2)
To ensure a robust evaluation and avoid training bias, we reserve orthologs of mouse chr1 and chr2 as our test set.

```bash
# Extract peaks mapped to mouse chr1 and chr2
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz \
  | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' \
  > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.narrowPeak
```

## 3. Filter Positive Test Set
We filter the Test peak set against our log-signal file to get a file with the OCR coordinates and respective logged signal values.
```bash
# Filter test set peaks using the log signal values file
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py \
  --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed \
  --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.narrowPeak \
  --peakNameCol 3 \
  --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos/rat_liver_TEST_500bp.bed
```

## 4. Identify Test 3 Peak Set
We identify regions that are accessible in Rat but not in Mouse. This is achieved by subtracting all reproducible and non-reproducible mouse peaks from our mapped rat peaks.
```bash
# Subtract mouse peaks from mapped rat peaks
# -A ensures that any overlap results in the removal of the entire rat peak
bedtools subtract -A \
  -a /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.narrowPeak \
  -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz \
  > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.narrowPeak
```

## 5. Filter Test 3 Set
We filter the Test 3 peak set against our log-signal file to get a file with the OCR coordinates and respective logged signal values.
```bash
# Filter test 3 set peaks using the log signal values file
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py \
  --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed \
  --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.narrowPeak \
  --peakNameCol 3 \
  --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/rat_liver_TEST_500bp.bed
```
