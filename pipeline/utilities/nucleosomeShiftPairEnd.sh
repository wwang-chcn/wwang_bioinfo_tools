#!/usr/bin/env bash

# USAGE: $0 <reads.bed>

function print_help {
	echo "USAGE: $0 <fragments.bed>"
}


awk '
{
	mid=($2+$3)/2
	if(mid-37>=0) {printf $1"\t%d\t%d\t"$4"\t"$5"\t"$6"\n", mid-37, mid+37}
	else printf $1"\t0\t%d\t"$4"\t"$5"\t"$6"\n", mid+37
}' $1 | sort -k1,1 -k2,2g > ${1::(${#1}-4)}_shift.bed

