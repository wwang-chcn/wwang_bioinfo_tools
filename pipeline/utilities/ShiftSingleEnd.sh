#! /bin/bash

# Nov-1-2018

# USAGE: $0 <reads.bed> <fragments_length>

function print_help {
        echo "USAGE: $0 <reads.bed>"
}


awk -v fl=${2} '
{if($6=="+") printf $1"\t%d\t%d\t"$4"\t"$5"\t+\n", $2+fl/4, $2+fl/4*3; else {if($3-fl/4>0) {if($3-fl/4*3<0) printf $1"\t0\t%d\t"$4"\t"$5"\t-\n", $3-fl/4; else printf $1"\t%d\t%d\t"$4"\t"$5"\t-\n", $3-fl/4*3, $3-fl/4}}}' $1 | sort -k1,1 -k2,2n | intersectBed -a - -b <(awk '{print $1"\t0\t"$2}' $3) -wa -f 1.00 > ${1::(${#1}-4)}_shift.bed
