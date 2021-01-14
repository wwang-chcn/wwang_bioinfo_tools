#!/bin/bash

# USAGE: $0 <mcall.out> <outputName.wig> [threshold=0]

function print_help {
	echo "USAGE: $0 <mcall.out> <outputName.wig> [threshold=0]"
}

if [[ -e $1 ]] && [[ -n $2 ]] 
then
	if [[ -n $3 ]]	# has threshold value
	then
		if [[ $3 =~ ^[1-9][0-9]?$ ]]	# input is a positive interge(right input)
		then
			threshold=$3
		else
			print_help; exit 2
		fi
	else
		threshold=1
	fi
else
	print_help; exit 1
fi

awk -v threshold=$threshold 'BEGIN{chr=""} NR>1{if($1!=chr) {chr=$1; printf "variableStep chrom=%s span=2\n", chr}; if($5>=threshold) {printf "%d\t%s\n", $2+1, $4}}' $1 > $2
awk -v threshold=$threshold 'BEGIN{chr=""} NR>1{if($1!=chr) {chr=$1; printf "variableStep chrom=%s span=2\n", chr}; if($5>=threshold) {printf "%d\t1\n", $2+1}}' $1 > ${2::(${#i}-4)}_control.wig
