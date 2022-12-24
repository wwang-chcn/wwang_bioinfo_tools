#! /bin/bash

function print_help {
    echo "$0 <input> <name>"
    echo "Generated cutsites Bed format files (${name}_cutsites.bed, ${name}_cutsites_plus.bed, and ${name}_cutsites_minus.bed) from fragments file (${input})."
}

function cut_sites {
    if [[ ! -e ${input} ]]; then
        print_help
        exit 1
    fi
    if ([ ! -e ${name}_cutsites.bed ] && [ ! -e ${name}_cutsites.bb ]) || ([ ! -e ${name}_cutsites_plus.bed ] && [ ! -e ${name}_cutsites_plus.bb ]) || ([ ! -e ${name}_cutsites_minus.bed ] && [ ! -e ${name}_cutsites_minus.bb ]); then
        awk '{print $1"\t"$2"\t"$2+1"\tcut_site\t0\t+\n"$1"\t"$3-1"\t"$3"\tcut_site\t0\t-"}' ${input} | sort -k1,1 -k2,2n > ${name}_cutsites.bed
        awk '{if($6=="+") print $0}' ${name}_cutsites.bed | sort -k1,1 -k2,2n > ${name}_cutsites_plus.bed
        awk '{if($6=="-") print $0}' ${name}_cutsites.bed | sort -k1,1 -k2,2n > ${name}_cutsites_minus.bed
    fi
}

input=${1}
name=${2}
cut_sites