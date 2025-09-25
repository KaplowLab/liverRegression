# MAP NON-MOUSE ALL PEAKS TO MOUSE

## RAT
sbatch -p pfen1 -w compute-1-11 -o ~/rat_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Rattus_norvegicus -t Mus_musculus -o /scratch/azstephe/liver/ratToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/rat_liver_pos_ALL.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

## MACAQUE
sbatch -p pfen1 -w compute-1-11 -o ~/mac_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Macaca_mulatta -t Mus_musculus -o /scratch/azstephe/liver/macaqueToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/macaque_liver_pos_ALL_genBankNames.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

## COW
sbatch -p pfen1 -w compute-1-11 -o ~/cow_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Bos_taurus -t Mus_musculus -o /scratch/azstephe/liver/cowToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

## PIG
sbatch -p pfen1 -w compute-1-11 -o ~/pig_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Sus_scrofa -t Mus_musculus -o /scratch/azstephe/liver/pigToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/pig_liver_pos_ALL_genBankNames.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

# CONVERT PIG TO GENBANKNAMES
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/raw/pig_liver_pos_ALL.narrowPeak --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigGenome/GenbankNameToChromNameSusScr3.txt --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/raw/pig_liver_pos_ALL_genBankNames.narrowPeak

# MAP MOUSE VAL TO THE 4 OTHER SPECIES
bash /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Mus_musculus -t Macaca_mulatta,Rattus_norvegicus,Bos_taurus,Sus_scrofa -o /scratch/azstephe/liver/mouseTo_VAL/ -b /home/azstephe/liverRegression/regression_liver/data/raw/mouse_liver_pos_VAL.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

# ################################################################################################################################################
# VAL SET CREATION
# /home/azstephe/liverRegression/regression_liver/data/splits/logPos/
# /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/
# /home/azstephe/liverRegression/regression_liver/data/splits/log_val2/
# /home/azstephe/liverRegression/regression_liver/data/splits/log_val3/
# ################################################################################################################################################

# VAL SET MOUSE COORDINATES, NON-MOUSE PEAKNAMES

## MACAQUE
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/macaque_liver_pos_ALL_genBankNames.Macaca_mulattaToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr8[[:space:]]|chr9[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_VAL.bed 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_VAL.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/logPos/macaque_liver_VAL.narrowPeak

## RAT
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr8[[:space:]]|chr9[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_VAL.bed 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_VAL.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/logPos/rat_liver_VAL.narrowPeak

## COW
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/cow_liver_pos_ALL_genBankNames.Bos_taurusToMus_musculus.HALPER.narrowPeak.gz  | grep -E '^(chr8[[:space:]]|chr9[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_VAL.bed 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_VAL.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/logPos/cow_liver_VAL.narrowPeak

## PIG
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/pig_liver_pos_ALL_genBankNames.Sus_scrofaToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr8[[:space:]]|chr9[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_VAL.bed 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_VAL.bed  --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/logPos/pig_liver_VAL.narrowPeak

# ################################################################################################################################################

# VAL 1

1. convert to [mouse val mapped to non-mouse] chromName if necessary
2. [mouse val mapped to non-mouse] subtract [non-mouse peaks rep & nonrep] > [val mouse peaks that are nonenhancer in non-mouse]
3. summit center extend to 500 bp *my way of extending may lead to negative indexes or out of sequence bounds*

## MACAQUE
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak.gz --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueGenome/GenbankNameToChromNameRheMac8.txt --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak_chromName.gz 

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak_chromName.gz -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-macs2_pooled/execution/Mac1_Lv_1_PFE116A29_S8_R1_001.merged.nodup.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/macaque_liver_VAL_raw.bed 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/macaque_liver_VAL_raw.bed > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/macaque_liver_VAL_500bp.bed 

## RAT
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToRattus_norvegicus.HALPER.narrowPeak.gz -b /projects/pfenninggroup/MPRA/Irene/rats/atac-pipeline-output/Liver/atac/a0da491a-cfc9-4701-97ab-c8cbfc0ef7bb/call-call_peak_pooled/execution/glob-2e6c87fc90e45fa5dbda88934d454cd3/R2-Liv_R1_001.trim.merged.nodup.no_chrM.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/rat_liver_VAL_raw.bed 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/rat_liver_VAL_raw.bed > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/rat_liver_VAL_500bp.bed 

## COW
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToBos_taurus.HALPER.narrowPeak.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/genomes/ChromNameToGenbankName_Btau_5.0.1 --chromNameDictReverse --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToBos_taurus.HALPER.narrowPeak_chromName.gz 

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToBos_taurus.HALPER.narrowPeak_chromName.gz -b /home/azstephe/liverRegression/regression_liver/data/raw/combinedCow_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/cow_liver_VAL_raw.bed 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/cow_liver_VAL_raw.bed > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/cow_liver_VAL_500bp.bed 

## PIG
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToSus_scrofa.HALPER.narrowPeak.gz --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigGenome/GenbankNameToChromNameSusScr3.txt --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToSus_scrofa.HALPER.narrowPeak_chromName.gz 

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToSus_scrofa.HALPER.narrowPeak_chromName.gz -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigDNase/Liver/atac/568838bc-4846-4491-8787-45b51e65bfee/peak/pooled-rep/ATAC_Liver_P348_R1.trim.merged.nodup.no_chrM.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/pig_liver_VAL_raw_chromName.bed 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/pig_liver_VAL_raw_chromName.bed > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/pig_liver_VAL_500bp.bed 

# ################################################################################################################################################

# VAL 2 

1. [non-mouse peaks mapped to mouse filtered for val set] intersect [raw mouse peaks] > [non-mouse val peaks, enhancer in mouse] 
2. [non-mouse val peaks, enhancer in mouse] filterPeakName with peaks from [set with transformed signal value] > [non-mouse val coords + peaks, orthologs enhancer in mouse] 

## MACAQUE
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_VAL.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseEnhancer_wawb_VAL.narrowPeak 

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseEnhancer_wawb_VAL.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val2/macaque_liver_VAL.narrowPeak 

## RAT
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_VAL.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb_VAL.narrowPeak 

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb_VAL.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val2/rat_liver_VAL.narrowPeak

## PIG
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_VAL.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseEnhancer_wawb_VAL.narrowPeak 

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseEnhancer_wawb_VAL.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val2/pig_liver_VAL.narrowPeak 

## COW
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_VAL.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseEnhancer_wawb_VAL.narrowPeak 

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseEnhancer_wawb_VAL.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val2/cow_liver_VAL.narrowPeak 

# ################################################################################################################################################

# VAL 3

1. [non-mouse all peaks mapped to mouse] subtract [raw mouse peaks (rep & nonrep)] > [non-mouse val peaks, nonenhancer in mouse]
2. [val set: non-mouse peaks that map to chr 8, 9 on mouse]
3. [non-mouse val peaks, nonenhancer in mouse] filterPeakName with peaks from [set with transformed signal value] > [non-mouse val coords + peaks, orthologs nonenhancer in mouse]

## MACAQUE

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/macaque_liver_pos_ALL_genBankNames.Macaca_mulattaToMus_musculus.HALPER.narrowPeak.gz -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseNon_ALL.bed

grep -E '^(chr8[[:space:]]|chr9[[:space:]])' macaqueToMouse_liver_macaqueEnhancer_mouseNon_ALL.bed > macaqueToMouse_liver_macaqueEnhancer_mouseNon_VAL3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseNon_VAL3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val3/macaque_liver_VAL.narrowPeak 

## RAT

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_ALL.bed

grep -E '^(chr8[[:space:]]|chr9[[:space:]])' ratToMouse_liver_ratEnhancer_mouseNon_ALL.bed > ratToMouse_liver_ratEnhancer_mouseNon_VAL3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/splits/ratMouse/ratToMouse_liver_ratEnhancer_mouseNon_val.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val3/rat_liver_VAL.narrowPeak 

## COW

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/cow_liver_pos_ALL_genBankNames.Bos_taurusToMus_musculus.HALPER.narrowPeak.gz -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseNon_ALL.bed

grep -E '^(chr8[[:space:]]|chr9[[:space:]])' cowToMouse_liver_cowEnhancer_mouseNon_ALL.bed > cowToMouse_liver_cowEnhancer_mouseNon_VAL3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/splits/cowMouse/cowToMouse_liver_cowEnhancer_mouseNon_val.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val3/cow_liver_VAL.narrowPeak

## PIG
bedtools subtract -A -a pig_liver_pos_ALL_genBankNames.Sus_scrofaToMus_musculus.HALPER.narrowPeak.gz -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseNon_ALL.bed

grep -E '^(chr8[[:space:]]|chr9[[:space:]])' pigToMouse_liver_pigEnhancer_mouseNon_ALL.bed > pigToMouse_liver_pigEnhancer_mouseNon_VAL3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/splits/pigMouse/pigToMouse_liver_pigEnhancer_mouseNon_val.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/log_val3/pig_liver_VAL.narrowPeak
 
# ################################################################################################################################################
# TEST SET CREATION
# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/
# /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/
# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/
# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/
# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/
# ################################################################################################################################################

# LOG TEST SET

1. get chr1, 2 from [non-mouse all mapped to mouse] > [non-mouse mapped to mouse TEST]
2. [non-mouse mapped to mouse TEST] filterPeakName with peaks from [set with transformed signal value] > [non-mouse test set]

## MACAQUE
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/macaque_liver_pos_ALL_genBankNames.Macaca_mulattaToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/macaque_liver_TEST_500bp.bed

## RAT
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/rat_liver_TEST_500bp.bed

## COW 
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/cow_liver_pos_ALL_genBankNames.Bos_taurusToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_TEST.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/cow_liver_TEST_500bp.bed

## PIG

zcat /home/azstephe/liverRegression/regression_liver/data/mapped/pig_liver_pos_ALL_genBankNames.Sus_scrofaToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_TEST.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/pig_liver_TEST_500bp.bed

# ################################################################################################################################################

# TEST NEG

## MACAQUE
awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rheMac8Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonMacaqueLiver_test_andRat_andCow_andPig.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/macaque_liver_TEST_500bp.bed

## RAT
awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rn6Fixed_summitExtendedMin50Max2XProtect5_nonRatLiver_test_andMacaque_andCow_andPig.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/rat_liver_TEST_500bp.bed 

## COW
awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_test_andMacaque_andRat_andPig.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/cow_liver_TEST_refSeq500bp.bed 

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/cow_liver_TEST_refSeq_500bp.bed --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/cow_liver_TEST_500bp.bed 

## PIG
awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_susScr3Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonPigLiver_test_andMacaque_andRat_andCow.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/pig_liver_TEST_500bp.bed 

# ################################################################################################################################################

# TEST 1

1. convert to [mouse test mapped to non-mouse] chromName if necessary
2. [mouse test mapped to non-mouse] subtract [non-mouse peaks rep & nonrep] > [test mouse peaks that are nonenhancer in non-mouse]
3. summit center extend to 500 bp *my way of extending may lead to negative indexes or out of sequence bounds*

## MACAQUE

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak.gz --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueGenome/GenbankNameToChromNameRheMac8.txt --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak_chromName.gz 

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak_chromName.gz -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-macs2_pooled/execution/Mac1_Lv_1_PFE116A29_S8_R1_001.merged.nodup.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/macaque_liver_TEST_raw.bed 

 awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/macaque_liver_TEST_raw.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/macaque_liver_TEST_500bp.bed 

## RAT
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToRattus_norvegicus.HALPER.narrowPeak.gz -b /projects/pfenninggroup/MPRA/Irene/rats/atac-pipeline-output/Liver/atac/a0da491a-cfc9-4701-97ab-c8cbfc0ef7bb/call-call_peak_pooled/execution/glob-2e6c87fc90e45fa5dbda88934d454cd3/R2-Liv_R1_001.trim.merged.nodup.no_chrM.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/rat_liver_TEST_raw.bed 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/rat_liver_TEST_raw.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/rat_liver_TEST_500bp.bed 

## COW
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToBos_taurus.HALPER.narrowPeak.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/genomes/ChromNameToGenbankName_Btau_5.0.1 --chromNameDictReverse --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToBos_taurus.HALPER.narrowPeak_chromName.gz 

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/pooled-rep/ATAC_Liver_M08_R1.trim.merged.nodup.no_NC_006853.1.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/genomes/ChromNameToRefSeqName_Btau_5.0.1 --chromNameDictReverse --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/raw/combinedCow_rep_nonRep.narrowPeak.gz 

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToBos_taurus.HALPER.narrowPeak_chromName.gz -b /home/azstephe/liverRegression/regression_liver/data/raw/combinedCow_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/cow_liver_TEST_raw.bed 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/cow_liver_TEST_raw.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/cow_liver_TEST_500bp.bed 

## PIG
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToSus_scrofa.HALPER.narrowPeak.gz --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigGenome/GenbankNameToChromNameSusScr3.txt --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToSus_scrofa.HALPER.narrowPeak_chromName.gz 

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo/mouse_liver_pos_TEST.Mus_musculusToSus_scrofa.HALPER.narrowPeak_chromName.gz -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigDNase/Liver/atac/568838bc-4846-4491-8787-45b51e65bfee/peak/pooled-rep/ATAC_Liver_P348_R1.trim.merged.nodup.no_chrM.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/pig_liver_TEST_raw_chromName.bed 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/pig_liver_TEST_raw_chromName.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/pig_liver_TEST_500bp.bed 

# ################################################################################################################################################

# TEST2

1. [non-mouse peaks mapped to mouse filtered for test set] intersect [raw mouse peaks] > [non-mouse test peaks, enhancer in mouse] 
2. [non-mouse test peaks, enhancer in mouse] filterPeakName with peaks from [set with transformed signal value] > [non-mouse test coords + peaks, orthologs enhancer in mouse]

## RAT
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb.narrowPeak 

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/rat_liver_TEST_500bp.bed

## MACAQUE
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseEnhancer_wawb.narrowPeak

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/macaque_liver_TEST_500bp.bed

## COW
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseEnhancer_wawb.narrowPeak

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/cow_liver_TEST_500bp.bed

## PIG
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseEnhancer_wawb.narrowPeak

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/pig_liver_TEST_500bp.bed

# ################################################################################################################################################

# TEST3

1. [non-mouse peaks mapped to mouse TEST] subtract [raw mouse peaks (rep & nonrep)] > [non-mouse test peaks, nonenhancer in mouse]
2. [non-mouse test peaks, nonenhancer in mouse] filterPeakName with peaks from [set with transformed signal value] > [non-mouse test coords + peaks, orthologs nonenhancer in mouse]

## RAT
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/rat_liver_TEST_500bp.bed

## MACAQUE
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseNon_TEST3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/macaque_liver_TEST_500bp.bed

## COW
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseNon_TEST3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/cow_liver_TEST_500bp.bed

## PIG
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseNon_TEST3.bed

### LOG
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/pig_liver_TEST_500bp.bed

# ################################################################################################################################################
# COW + PIG CREATION
# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test4/
# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test5/
# /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/ (val + train)
# ################################################################################################################################################

### convert to chromName [train, val, test sets positives]
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_inPigLiverLoose_nonMouseLiver_nonRatLiver_nonMacaqueLiver_train.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_TRAIN_chromName.bed.gz

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_inPigLiverLoose_nonMouseLiver_nonRatLiver_nonMacaqueLiver_valid.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_VAL_chromName.bed.gz

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_inPigLiverLoose_nonMouseLiver_nonRatLiver_nonMacaqueLiver_test.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --chromNameDictReverse --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test4/cow_pig_liver_pos_mouse_macaque_rat_closed_TEST_chromName.bed.gz 

### convert to chromName [train, val, test sets negatives]
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllNonLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_inRatLiverLoose_inMacaqueLiverLoose_nonPigLiver_train.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_TRAIN_chromName.bed.gz

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllNonLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_inRatLiverLoose_inMacaqueLiverLoose_nonPigLiver_valid.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_VAL_chromName.bed.gz

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllNonLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_inRatLiverLoose_inMacaqueLiverLoose_nonPigLiver_test.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --chromNameDictReverse --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test5/cow_pig_liver_neg_mouse_macaque_rat_open_TEST_chromName.bed.gz 

# get cow log positive values
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_TRAIN_chromName.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_TRAIN_chromName_500bp.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_VAL_chromName.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_VAL_chromName_500bp.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName  /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName  /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test4/cow_pig_liver_pos_mouse_macaque_rat_closed_TEST_chromName.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test4/cow_pig_liver_pos_TEST_500bp.bed 

# summit center and extend negatives to 500 bp
awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_TRAIN_chromName.bed > /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_TRAIN_chromName_500bp.bed

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_VAL_chromName.bed > /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_VAL_chromName_500bp.bed

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test5/cow_pig_liver_neg_mouse_macaque_rat_open_TEST_chromName.bed > /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test5/cow_pig_liver_neg_TEST_500bp.bed 

# ################################################################################################################################################
# 2KB AND EQN SET CREATION
# ################################################################################################################################################

1. Basically filters /home/azstephe/liverRegression/regression_liver/data/ladder_qn/ instead of /home/azstephe/liverRegression/regression_liver/data/log/

# Run:
bash /home/azstephe/liverRegression/regression_liver/scripts/make_2kb.sh
bash /home/azstephe/liverRegression/regression_liver/scripts/make_eqn.sh

# ################################################################################################################################################
# USED THESE SETS FOR EVAL
# ################################################################################################################################################

# /home/azstephe/liverRegression/regression_liver/data/splits/logPos/*VAL* 

# /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/ 

# /home/azstephe/liverRegression/regression_liver/data/splits/log_val2/ 

# /home/azstephe/liverRegression/regression_liver/data/splits/log_val3/ 

# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/ 

# /home/azstephe/liverRegression/regression_liver/data/test_splits/neg/ 

# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test1/ 

# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test2/ 

# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test3/ 

# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test4/cow_pig_liver_pos_TEST_500bp.bed 

# /home/azstephe/liverRegression/regression_liver/data/test_splits/log_test5/cow_pig_liver_neg_TEST_500bp.bed 

# /home/azstephe/liverRegression/regression_liver/data/splits/logPos/*TRAINONLY.narrowPeak* 
