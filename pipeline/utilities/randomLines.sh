#! /bin/bash

n=`wc -l ${1} | cut -f 1 -d " "`
c=`bc -l <<< "${2} / ${n}"`
awk -v ratio=$c 'BEGIN{srand(6666)} {if(rand() < ratio) print $0}' ${1} > ${3}