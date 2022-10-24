#! /bin/bash

# ----- misc funcs -----

function print_help {
    echo "$0 <name> <genomeVersion> <normalizationFactor>"
    echo "<normalizationFactor> is control the scale factor parameter for genomeCoverageBed step."
    echo "    The value could be 'true', 'false', or a positive integer."
    echo "    A integer value means how many reads should be sampling to."
    echo "    'true' means sampling to 1M reads."
    echo "    'false' means do not perform sampling."
}

function is_positive_integer {
    if [[ ${1} =~ ^[\-0-9]+$ ]] && (( ${1} > 0)); then
        func_res="true"
    else
        func_res="false"
    fi
}

function normalizationFactor_handling {
    case ${1} in 
        true ) normalizationFactor=1000000;;
        false ) normalizationFactor=false;;
        * ) normalizationFactor=${1} && is_positive_integer ${normalizationFactor} ;;
    esac

    case ${func_res} in
        true ) ;;
        false ) echo "Invalid normalizationFactor parameters:" ${normalizationFactor} && print_help;;
    esac
}

# ----- parameters -----

if [[ $# -lt 4 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

name=${1}
genomeVersion=${2}

normalizationFactor_handling ${3}

# get MY_PATH (/wwang_bioinfo_tools/pipeline/utilities)
MY_PATH="`dirname \"$0\"`"

#----- pileup -----

if [[ ! -e 2_signal/${name}.bw ]]; then
    cd 2_signal && \
    awk '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$6}' ${name}_raw_reads.bed | sort -S 1% -k1,1 -k2,2n | uniq | awk '{print $1"\t"$2"\t"$3"\tReads"NR"\t"$5"\t"$6}' > ${name}_reads.bed
    ${MY_PATH}/nucleosomeShiftSingleEnd.sh ${name}_reads.bed && \
    n=`wc -l ${name}_reads_shift.bed | cut -f 1 -d " "` && \
    c=`bc -l <<< "1000000 / $n"` && \
    genomeCoverageBed -bga -scale $c -i ${name}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name}_reads_shift.bdg && \
    ${MY_PATH}/bdg2bw.sh ${name}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
    rm ${name}_reads_shift.bed ${name}_reads_shift.bdg
    cd ..
fi