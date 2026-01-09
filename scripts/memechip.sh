#!/bin/bash

# Description: Submits MEME-ChIP as an SBATCH job
# Usage bash memechip.sh <output director> <foreground fasta> <background fasta>

outdir=$1
infasta=$2
backgroundfasta=$3

~/anaconda2/bin/meme-chip -o $outdir -db ~/MotifData/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme -neg $backgroundfasta $infasta

