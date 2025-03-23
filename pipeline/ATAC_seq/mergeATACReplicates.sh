#! /bin/bash

# ----- functions -----

function print_help {
    echo "USAGE: $0 <name> <genomeVersion> <replicates,+> [downSampling]"
    echo "name: the sample name"
    echo "genomeVersion: the genome version processed for these samples"
    echo "replicates: comma separated list of replicates"
    echo "downSampling: target sampling fragments for each sample. optional"
}

function bedToBigWig {
    name=${1}
    genomeVersion=${2}
    n=`wc -l ${name}_OCR_fragments.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"`
    genomeCoverageBed -bga -scale $c -i ${name}_OCR_fragments.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_OCR_fragments.bdg
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_OCR_fragments.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}
    rm ${name}_OCR_fragments.bdg
}

function compress_bed {
    bedFile=${1}
    genomeVersion=${2}
    col=`head -1 ${bedFile} | awk '{print NF}'`
    plus=`bc <<< "$col -3"`
    intersectBed -a ${bedFile} -b <(awk '{print $1"\t0\t"$2}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes) -wa -f 1.00 | sort -k1,1 -k2,2n > ${bedFile}.tmp
    bedToBigBed -type=bed3+${plus} ${bedFile}.tmp ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes ${bedFile::(${#bedFile}-2)}b
    rm ${bedFile} ${bedFile}.tmp
}

function peak_calling {
    name=${1}
    genomeVersion=${2}
    
    mkdir -p ../4_peaks
    cd 4_peaks
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes`
    if [[ ! -e ../3_merged_signal/${name}_OCR_fragments.bed ]]; then
        bigBedToBed ../3_merged_signal/${name}_OCR_fragments.bb ../3_merged_signal/${name}_OCR_fragments.bed
    fi
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
    macs3 callpeak -f BEDPE -t ../3_merged_signal/${name}_OCR_fragments.bed -n ${name} -g ${chromsize} --keep-dup all 2>&1 >>/dev/null | tee MACS_${name}.out
    cd ..
}

function clearning_up {
    for i in $@; do
        if [[ -e 2_signal/OCR/${i}_uniq_OCR_fragments.bb ]]; then
            rm 2_signal/OCR/${i}_uniq_OCR_fragments.bed
        fi
    done
    if [[ -e ../3_merged_signal/${name}_OCR_fragments.bed ]]; then
        rm ../3_merged_signal/${name}_OCR_fragments.bed
    fi
}

if [[ $# -lt 3 ]]; then
    echo "No enought parameters!"
    print_help
    exit 1
fi

if [[ $# -gt 4 ]]; then
    echo "Too many parameters!"
    print_help
    exit 1
fi

MY_PATH="`readlink -f $(dirname \"$0\")`"

mkdir -p 3_merged_signal/ 4_peaks/


name=${1}
genomeVersion=${2}
IFS=',' read -r -a reads_file_array <<< ${3}
for i in ${reads_file_array[@]}; do
    if [[ ! -e 2_signal/OCR/${i}_OCR_fragments.bed ]]; then
        if [[ ! -e 2_signal/OCR/${i}_OCR_fragments.bb ]]; then
            echo "fragments for sample ${i} do not exist, exit!"
            exit 1
        fi
        bigBedToBed 2_signal/OCR/${i}_OCR_fragments.bb 2_signal/OCR/${i}_OCR_fragments.bed
    fi
    fragments_array+=("2_signal/OCR/${i}_OCR_fragments.bed")
done

cat ${fragments_array[@]} | sort -k1,1 -k2,2n > 3_merged_signal/${name}_OCR_fragments.bed.tmp

cd 3_merged_signal/
if [[ $# -lt 4 ]]; then
    mv ${name}_OCR_fragments.bed.tmp ${name}_OCR_fragments.bed
else
    n=`wc -l ${name}_OCR_fragments.bed.tmp | cut -f 1 -d " "`
    ratio=`bc -l <<< "${4} / ${n}"`
    awk -v ratio=${ratio} 'BEGIN{srand(1006)} {if(rand()<ratio) print}' ${name}_OCR_fragments.bed.tmp > ${name}_OCR_fragments.bed
    rm ${name}_OCR_fragments.bed.tmp
fi
bedToBigWig ${name} ${genomeVersion}
compress_bed ${name}_OCR_fragments.bed ${genomeVersion}
cd ..

peak_calling ${name} ${genomeVersion}
clearning_up
