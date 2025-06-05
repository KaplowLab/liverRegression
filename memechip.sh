#!/bin/bash

outdir=$1
infasta=$2
backgroundfasta=$3

/home/ikaplow/anaconda2/bin/meme-chip -o $outdir -db /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme -neg $backgroundfasta $infasta

#/home/ikaplow/anaconda2/bin/meme-chip -o $outdir -db /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme $infasta
