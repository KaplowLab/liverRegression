################## RAT
# map all rat liver peaks to mouse
sbatch -p pfen1 -w compute-1-11 -o ~/rat_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Rattus_norvegicus -t Mus_musculus -o /scratch/azstephe/liver/ratToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/rat_liver_pos_ALL.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

>Mapped 17477 peaks out of 21531 total peaks

# TEST3
        # grep chr1, 2
zcat /home/azstephe/liverRegression/regression_liver/data/mapped/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed

# TEST
python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/rat_liver_TEST_500bp.bed

        # rat test set mapped to mouse subtract mouse rep+nonrep

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.bed

# rat enhancers, mouse non ALL chrom (including val set)
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_ALL.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test3/rat_liver_TEST_500bp.bed

#VAL3

grep -E '^(chr8[[:space:]]|chr9[[:space:]])' ratToMouse_liver_ratEnhancer_mouseNon_ALL.bed > ratToMouse_liver_ratEnhancer_mouseNon_VAL3.bed

# TEST2

bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb.narrowPeak 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test2/rat_liver_TEST_500bp.bed

zcat /home/azstephe/liverRegression/regression_liver/data/mapped/rat_liver_pos_ALL.Rattus_norvegicusToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr8[[:space:]]|chr9[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_VAL.bed

############### MACAQUE
# map all macaque liver peaks to mouse
sbatch -p pfen1 -w compute-1-11 -o ~/mac_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Macaca_mulatta -t Mus_musculus -o /scratch/azstephe/liver/macaqueToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/macaque_liver_pos_ALL_genBankNames.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

>Mapped 16486 peaks out of 32697 total peaks.

# TEST

zcat /home/azstephe/liverRegression/regression_liver/data/mapped/macaque_liver_pos_ALL_genBankNames.Macaca_mulattaToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/macaque_liver_TEST_500bp.bed

# TEST3

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseNon_TEST3.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test3/macaque_liver_TEST_500bp.bed

# TEST2

bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseEnhancer_wawb.narrowPeak

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test2/macaque_liver_TEST_500bp.bed

############### COW
# map all cow liver peaks to mouse
        # get chr13 peaks
grep CM000189.7 /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_chr13.narrowPeak
        # isolate nonchr13 peaks
awk '$1!="CM000189.7"' /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_nochr13.narrowPeak
        # split on last instance of chr12 - insertion location for chr13
awk '$1=="CM000188.7"{seen=1; last=NR} {lines[NR]=$0} END{for(i=1;i<=NR;i++){if(i<=last) print lines[i] > "/home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_nochr13.narrowPeak_part1"; else print lines[i] > "/home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_nochr13.narrowPeak_part2"}}' /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_nochr13.narrowPeak
        # concat all 3 files
cat /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_nochr13.narrowPeak_part1 /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_chr13.narrowPeak /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames_nochr13.narrowPeak_part2 > /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames.narrowPeak

sbatch -p pfen1 -w compute-1-11 -o ~/cow_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Bos_taurus -t Mus_musculus -o /scratch/azstephe/liver/cowToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/cow_liver_pos_ALL_genBankNames.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

# all cow open, mouse closed peaks
bedtools subtract -A -a cow_liver_pos_ALL_genBankNames.Bos_taurusToMus_musculus.HALPER.narrowPeak.gz -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseNon_ALL.bed

# val3
grep -E '^(chr8[[:space:]]|chr9[[:space:]])' cowToMouse_liver_cowEnhancer_mouseNon_ALL.bed > cowToMouse_liver_cowEnhancer_mouseNon_VAL3.bed

# ENTIRE TEST SET
zcat cow_liver_pos_ALL_genBankNames.Bos_taurusToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > cowToMouse_liver_cowEnhancer_TEST.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/cow_liver_TEST_500bp.bed
        # TEST3
bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseNon_TEST3.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test3/cow_liver_TEST_500bp.bed

        # TEST2
bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseEnhancer_wawb.narrowPeak

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test2/cow_liver_TEST_500bp.bed

        # TEST1
bash /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Mus_musculus -t Bos_taurus -o /scratch/azstephe/liver/mouseToCow/ -b /home/azstephe/liverRegression/regression_liver/data/raw/mouse_liver_pos_TEST.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

####################### PIG

# convert pig to genBankNames
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/raw/pig_liver_pos_ALL.narrowPeak --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigGenome/GenbankNameToChromNameSusScr3.txt --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/raw/pig_liver_pos_ALL_genBankNames.narrowPeak


# map all pig liver peaks to mouse
sbatch -p pfen1 -w compute-1-11 -o ~/pig_mouse.o /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Sus_scrofa -t Mus_musculus -o /scratch/azstephe/liver/pigToMouse/ -b /home/azstephe/liverRegression/regression_liver/data/raw/pig_liver_pos_ALL_genBankNames.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

# all pig open, mouse closed peaks
bedtools subtract -A -a pig_liver_pos_ALL_genBankNames.Sus_scrofaToMus_musculus.HALPER.narrowPeak.gz -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseNon_ALL.bed

#pig VAL3
grep -E '^(chr8[[:space:]]|chr9[[:space:]])' pigToMouse_liver_pigEnhancer_mouseNon_ALL.bed > pigToMouse_liver_pigEnhancer_mouseNon_VAL3.bed

zcat pig_liver_pos_ALL_genBankNames.Sus_scrofaToMus_musculus.HALPER.narrowPeak.gz | grep -E '^(chr1[[:space:]]|chr2[[:space:]])' > pigToMouse_liver_pigEnhancer_TEST.bed

bedtools intersect -wa -wb -a /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/mouse_liver_pos_ALL.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseEnhancer_wawb.narrowPeak

#TEST

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseEnhancer_wawb.narrowPeak --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test2/pig_liver_TEST_500bp.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_TEST.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_pos_LiuAll/pig_liver_TEST_500bp.bed

# TEST3

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_TEST.bed -b /home/azstephe/regression_liver/data/raw/combinedMouse_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseNon_TEST3.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/pig_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/pigToMouse_liver_pigEnhancer_mouseNon_TEST3.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/test_splits/log_LiuAll_test3/pig_liver_TEST_500bp.bed



##############################################################
VAL SET CREATION 9/10 ########################################

# map mouse val set to everything
bash /home/azstephe/liverRegression/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh -s Mus_musculus -t Macaca_mulatta,Rattus_norvegicus,Bos_taurus,Sus_scrofa -o /scratch/azstephe/liver/mouseTo_VAL/ -b /home/azstephe/liverRegression/regression_liver/data/raw/mouse_liver_pos_VAL.narrowPeak --halPath /home/ikaplow/RegulatoryElementEvolutionProject/src/hal/bin/halLiftover

# MACAQUE VAL1

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak.gz --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueGenome/GenbankNameToChromNameRheMac8.txt --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak_chromName.gz

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToMacaca_mulatta.HALPER.narrowPeak_chromName.gz -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-macs2_pooled/execution/Mac1_Lv_1_PFE116A29_S8_R1_001.merged.nodup.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/macaque_liver_VAL_raw.bed

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/macaque_liver_VAL_raw.bed > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/macaque_liver_VAL_500bp.bed

# MACAQUE VAL

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/macaque_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/macaqueToMouse_liver_macaqueEnhancer_VAL.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/val_splits/val/macaque_liver_VAL_500bp.bed

# RAT VAL1

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToRattus_norvegicus.HALPER.narrowPeak.gz -b /projects/pfenninggroup/MPRA/Irene/rats/atac-pipeline-output/Liver/atac/a0da491a-cfc9-4701-97ab-c8cbfc0ef7bb/call-call_peak_pooled/execution/glob-2e6c87fc90e45fa5dbda88934d454cd3/R2-Liv_R1_001.trim.merged.nodup.no_chrM.tn5.pooled.pval0.01.300K.bfilt.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/rat_liver_VAL_raw.bed

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/rat_liver_VAL_raw.bed > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/rat_liver_VAL_500bp.bed

# RAT VAL

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/rat_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/ratToMouse_liver_ratEnhancer_VAL.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/val_splits/val/rat_liver_VAL_500bp.bed

# COW VAL1

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToBos_taurus.HALPER.narrowPeak.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/genomes/ChromNameToGenbankName_Btau_5.0.1 --chromNameDictReverse --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToBos_taurus.HALPER.narrowPeak_chromName.gz

bedtools subtract -A -a /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToBos_taurus.HALPER.narrowPeak_chromName.gz -b /home/azstephe/liverRegression/regression_liver/data/raw/combinedCow_rep_nonRep.narrowPeak.gz > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/cow_liver_VAL_raw.bed

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/cow_liver_VAL_raw.bed > /home/azstephe/liverRegression/regression_liver/data/val_splits/val1/cow_liver_VAL_500bp.bed

# COW VAL

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mapped/cowToMouse_liver_cowEnhancer_VAL.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/val_splits/val/cow_liver_VAL_500bp.bed

# PIG VAL1

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToSus_scrofa.HALPER.narrowPeak.gz --chromNameDictFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigGenome/GenbankNameToChromNameSusScr3.txt --gzip --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouseTo_VAL/mouse_liver_pos_VAL.Mus_musculusToSus_scrofa.HALPER.narrowPeak_chromName.gz

# COW + PIG

        # convert to chromName train and val sets positives

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_inPigLiverLoose_nonMouseLiver_nonRatLiver_nonMacaqueLiver_valid.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_VAL_chromName.bed.gz

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_inPigLiverLoose_nonMouseLiver_nonRatLiver_nonMacaqueLiver_train.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_TRAIN_chromName.bed.gz

        # negs train and val
python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllNonLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_inRatLiverLoose_inMacaqueLiverLoose_nonPigLiver_valid.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_VAL_chromName.bed.gz

python /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/convertChromNames.py --bedFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllNonLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_inRatLiverLoose_inMacaqueLiverLoose_nonPigLiver_train.bed.gz --chromNameDictFileName /home/azstephe/liverRegression/regression_liver/data/ChromNameToRefSeqNameBtau_5.0.1.txt --gzip --chromNameDictReverse --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_TRAIN_chromName.bed.gz

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_VAL_chromName.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_VAL_chromName_500bp.bed

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_VAL_chromName.bed > /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_VAL_chromName_500bp.bed

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName /home/azstephe/liverRegression/regression_liver/data/log/cow_liver_pos_ALL.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_TRAIN_chromName.bed --peakNameCol 3 --outputFileName /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_pos_mouse_macaque_rat_closed_TRAIN_chromName_500bp.bed

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_TRAIN_chromName.bed > /home/azstephe/liverRegression/regression_liver/data/splits/cow_pig/cow_pig_liver_neg_mouse_macaque_rat_open_TRAIN_chromName_500bp.bed



