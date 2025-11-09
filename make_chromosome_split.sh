#/bin/bash

infile=$1
outfile=$2

awk -F"\t" '/chr1\t/ || /chr2\t/ {print $0}' "$infile"| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$7}' > "${outfile}_TEST.narrowPeak"

awk -F"\t" '/chr8\t/ || /chr10\t/ {print $0}' "$infile"| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$7}' > "${outfile}_VAL.narrowPeak"

awk -F"\t" '!/chr8\t/ || !/chr10\t/ {print $0}' "$infile" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$7}' >  "${outfile}_TRAIN.narrowPeak"

sed -i '/chr8\t/d' "${outfile}_TRAIN.narrowPeak"

sed -i '/chr10\t/d' "${outfile}_TRAIN.narrowPeak"

sed -i '/chr1\t/d' "${outfile}_TRAIN.narrowPeak"

sed -i '/chr2\t/d' "${outfile}_TRAIN.narrowPeak"
