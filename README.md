# liverRegression

## About

This repo contains the code ran for the project "Challenges in Predicting Chromatin Accessibility Differences between Species".
This includes a program for quantile normalization (quantile_normalize.py) and extended quantile normalization (EQN.py).


## Preprocessing the files

Start with the ATAC seq file which is organized with tab separated with columns 
@ gemini fill this in

Summit center the file 
`awk 'BEGIN{OFS="\t"} {print $1, $2 + $10 - 250, $2 + $10 + 250, $4, $7}' rat_liver_pos_ALL_notSC.bed > rat_liver_pos_ALL.bed`


### Test Set Creation

#### Set 1
Open in training species; closed in evaluation species

#### Set 2
Open in training species; open in evaluation species

#### Set 3
Closed in training species; open in evaluation species

