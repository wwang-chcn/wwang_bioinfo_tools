#!/usr/bin/env bash

# USAGE: $0 bdg gsize name

bedtools intersect -a $1 -b <(awk '{print $1"\t0\t"$2}' $2) -wa -f 1.00 | sort -S 5% -k1,1 -k2,2n > ${1}.tmp
bedGraphToBigWig ${1}.tmp $2 $3.bw
rm ${1}.tmp
