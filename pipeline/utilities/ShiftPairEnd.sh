#! /bin/bash

# Nov-1-2018

# USAGE: $0 <fragments.bed> <fragments_length> <chrom.sizes>

function print_help {
        echo "USAGE: $0 <fragments.bed> <fragments_length> <chrom.sizes>"
}


awk -v extsize=${2} '
{
        mid=($2+$3)/2
        printf $1"\t%d\t%d\n", mid-extsize/4, mid+extsize/4
}' $1 | awk '{if($2>=0) print; else print $1"\t0\t"$3}' | sort -k1,1 -k2,2n | intersectBed -a - -b <(awk '{print $1"\t0\t"$2}' $3) -wa -f 1.00 > ${1::(${#1}-4)}_shift.bed