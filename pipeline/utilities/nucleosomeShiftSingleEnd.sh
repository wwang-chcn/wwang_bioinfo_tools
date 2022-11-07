#!/usr/bin/env bash

# USAGE: $0 <reads.bed>

function print_help {
	echo "USAGE: $0 <reads.bed>"
}

shiftSize=37
extSize=73
endShiftSize=$((${shiftSize}+${extSize}))


awk -v shiftSize=$shiftSize -v endShiftSize=$endShiftSize '
{
	if($6=="+")
		printf $1"\t"$2+shiftSize"\t"$2+endShiftSize"\t"$4"\t"$5"\t"$6"\n"
	else
		{
			if($3-endShiftSize>0)
				printf $1"\t"$3-endShiftSize"\t"$3-shiftSize"\t"$4"\t"$5"\t"$6"\n"
			else
			    if($3-shiftSize>0) printf $1"\t0\t"$3-shiftSize"\t"$4"\t"$5"\t"$6"\n"
		}
}' $1 | sort -k1,1 -k2,2n > ${1::(${#1}-4)}_shift.bed

