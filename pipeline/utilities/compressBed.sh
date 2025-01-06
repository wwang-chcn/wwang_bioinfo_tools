#!/usr/bin/env bash

bedFile=${1}
genomeVersion=${2}

if [[ ! -e ${bedFile::(${#bedFile}-2)}b ]]; then
    col=`head -1 ${bedFile} | awk '{print NF}'` && \
    plus=`bc <<< "$col -3"` && \
    intersectBed -a ${bedFile} -b <(awk '{print $1"\t0\t"$2}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes) -wa -f 1.00 | sort -k1,1 -k2,2n > ${bedFile}.tmp && \
    bedToBigBed -type=bed3+${plus} ${bedFile}.tmp ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes ${bedFile::(${#bedFile}-2)}b && \
    rm ${bedFile}.tmp # && \
    rm ${bedFile}
fi