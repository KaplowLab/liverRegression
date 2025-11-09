#MAKE NARROWPEAK FILE FROM PEAKLIST 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName  /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_more_than_03_peaklist.csv --peakNameCol 0 --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_more_than_03.narrowPeak  

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName  /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_less_than_03_peaklist.csv --peakNameCol 0 --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_less_than_03.narrowPeak 

#GET FASTA FROM BED/NARROWPEAK FILE 

bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/mm10.fa -bed /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_more_than_03.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_more_than_03_sequences.fa 

bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/mm10.fa -bed /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_less_than_03.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_less_than_03_sequences.fa 

#RUN MEMECHIP 

/home/ikaplow/anaconda2/bin/meme-chip -o /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouseLess03Memechip -db /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme -neg /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_more_than_03_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_less_than_03_sequences.fa 

/home/ikaplow/anaconda2/bin/meme-chip -o /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouseLess03Memechip -db /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_less_than_03_sequences.fa 

/home/ikaplow/anaconda2/bin/meme-chip -o /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouseLess03Memechip -db /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouse_std_more_than_03_sequences.fa 

 

##MOUSE 

#  6459 mouse_between_1p4_2_peaklist.csv 

#  3785 mouse_greater_than_2p5_peaklist.csv 

#  2701 mouse_less_than_1p4_peaklist.csv 

##STRONG 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName  /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_peaklist.csv --peakNameCol 0 --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5.narrowPeak 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' mouse_greater_than_2p5.narrowPeak > mouse_greater_than_2p5_500bp.narrowPeak 

bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/mm10.fa -bed /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_500bp.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa 

##MEDIUM 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName  /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_peaklist.csv --peakNameCol 0 --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2.narrowPeak 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' mouse_between_1p4_2.narrowPeak > mouse_between_1p4_2_500bp.narrowPeak 

bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/mm10.fa -bed /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_500bp.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_sequences.fa 

##WEAK 

python ~/liverRegression/repos/OCROrthologPrediction/src/filterPeakName.py --unfilteredPeakFileName  /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort.bed --peakListFileName /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_peaklist.csv --peakNameCol 0 --outputFileName /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4.narrowPeak 

awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' mouse_less_than_1p4.narrowPeak > mouse_less_than_1p4_500bp.narrowPeak 

bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/mm10.fa -bed /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_500bp.narrowPeak > /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_sequences.fa 

##Strong vs weak peaks 

/home/ikaplow/anaconda2/bin/meme-chip -o /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseStrongVsWeakMemechip -db /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme -neg /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa 

(4797019) 

sbatch --mem=4G -o str_weak.o -J str_weak -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseStrongVsWeakMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_sequences.fa  

sbatch --mem=4G -o weak_str.o -J weak_str -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseWeakVsStrongMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa 

sbatch --mem=4G -o str_weak.o -J str_weak -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseStrongVsWeakMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_sequences.fa 

##Strong vs medium peaks 

/home/ikaplow/anaconda2/bin/meme-chip -o /home/azstephe/liverRegression/regression_liver/data/mouse_certainty/mouseLess03Memechip -db /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme -neg /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa 

 (4797018) 

sbatch --mem=4G -o str_med.o -J str_med -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseStrongVsMediumMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_sequences.fa 

sbatch --mem=4G -o med_str.o -J med_str -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseMediumVsStrongMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_sequences.fa /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa 
  
##Comparing weak, medium, strong peak motif enrichment 

 

sbatch --mem=4G -o med.o -J med -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseMediumMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_sequences.fa  

 

sbatch --mem=4G -o strong.o -J strong -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseStrongMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa 

 

sbatch --mem=4G -o weak.o -J weak -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseWeakMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_sequences.fa   

 

sbatch --mem=4G -o med.o -J med -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseMediumMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_between_1p4_2_sequences.fa  

 

sbatch --mem=4G -o strong.o -J strong -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseStrongMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_greater_than_2p5_sequences.fa 

 

sbatch --mem=4G -o weak.o -J weak -p gpu -n 1 --gres gpu:1 /home/azstephe/liverRegression/regression_liver/scripts/memechip.sh /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouseWeakMemechip /home/azstephe/liverRegression/regression_liver/data/mouse_strength/mouse_less_than_1p4_sequences.fa  
