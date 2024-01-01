#! /bin/bash

# USAGE: $0 <fragments.bed> <chrom.sizes> <name>

function print_help {
        echo "USAGE: $0 <fragments.bed> <chrom.sizes>  <name>"
}

awk '{
    mid=($2+$3)/2
    fraglen=$3-$2
    printf $1"\t%d\t%d\t%s\n", mid-fraglen/2, mid+fraglen/2, $4
}' $1 | awk '{if($2>=0) print; else print $1"\t0\t"$3"\t"$4}' | sort -k1,1 -k2,2n | intersectBed -a - -b <(awk '{print $1"\t0\t"$2}' $2) -wa -f 1.00 > ${1::(${#1}-4)}_shift.bed && \
n=`wc -l ${1::(${#1}-4)}_shift.bed | cut -f 1 -d " "` && \
c=`bc -l <<< "1000000 / $n"` && \
genomeCoverageBed -bga -scale $c -i ${1::(${#1}-4)}_shift.bed -g $2 | awk '{if($3>$2) print $0}' > ${1::(${#1}-4)}.bdg && \
bedtools intersect -a ${1::(${#1}-4)}.bdg -b <(awk '{print $1"\t0\t"$2}' $2) -wa -f 1.00 | sort -S 5% -k1,1 -k2,2n > ${1}.tmp
bedGraphToBigWig ${1}.tmp ${2} ${3}.bw
rm ${1::(${#1}-4)}_shift.bed ${1}.tmp