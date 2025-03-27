#! /bin/bash

# May-28-2021

function print_help {
    echo "USAGE: $0 <name> <genomeVersion> <replicates,+> [downSampling]"
    echo "name: the sample name"
    echo "genomeVersion: the genome version processed for these samples"
    echo "replicates: comma separated list of replicates"
    echo "downSampling: target sampling fragments for each sample. optional"
}


function bedToBigWig {
    echo "bedToBigWig ${1}" # debug
    fragment_length=`awk 'BEGIN{s=0} {s+=$3-$2} END{printf "%f", s/NR}' ${1}_fragments.bed`
    ${MY_PATH}/../utilities/ShiftPairEnd.sh ${1}_fragments.bed ${fragment_length} ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes
    n=`wc -l ${1}_fragments_shift.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"`
    genomeCoverageBed -bga -scale $c -i ${1}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${1}_fragments_shift.bdg
    ${MY_PATH}/../utilities/bdg2bw.sh ${1}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${1}
    rm ${1}_fragments_shift.bdg ${1}_fragments_shift.bed
}

function compress_bed {
    bedFile=${1}
    genomeVersion=${2}
    echo "compress_bed ${bedFile} ${genomeVersion}" # debug
    col=`head -1 ${bedFile} | awk '{print NF}'`
    plus=`bc <<< "$col -3"`
    intersectBed -a ${bedFile} -b <(awk '{print $1"\t0\t"$2}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes) -wa -f 1.00 | sort -k1,1 -k2,2n > ${bedFile}.compress.tmp
    bedToBigBed -type=bed3+${plus} ${bedFile}.compress.tmp ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes ${bedFile::(${#bedFile}-2)}b
    rm ${bedFile}.compress.tmp
}
function peak_calling {
    pc_name=${1}
    genomeVersion=${2}
    echo "peak_calling ${pc_name} ${genomeVersion}" # debug
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes`
    if [[ ! -e ${pc_name}_fragments.bed ]]; then
        bigBedToBed ${pc_name}_fragments.bb ${pc_name}_fragments.bed
    fi
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
    macs3 callpeak -f BEDPE -t ${pc_name}_fragments.bed -n ${pc_name} -g ${chromsize} --keep-dup all 2>&1 >>/dev/null | tee MACS_${pc_name}.out
}

if [[ $# -lt 4 ]]; then
    echo "No enought parameters!"
    print_help
    exit 1
fi

MY_PATH="`readlink -f $(dirname \"$0\")`"

mkdir -p 5_merged_sample/

name=${1}
genomeVersion=${2}

# ----- OCR fragments -----
# get fragment files
IFS=',' read -r -a reads_file_array <<< ${3}
for i in ${reads_file_array[@]}; do
    if [[ ! -e 2_signal/${i}_OCR_fragments.bed ]]; then
        if [[ ! -e 2_signal/${i}_OCR_fragments.bb ]]; then
            echo "fragments file for sample ${i} do not exist, exit!"
            exit 1
        fi
        bigBedToBed 2_signal/${i}_OCR_fragments.bb 2_signal/${i}_OCR_fragments.bed
    fi
    fragments_array+=("2_signal/${i}_OCR_fragments.bed")
done

# merge fragment files
cat ${fragments_array[@]} | sort -k1,1 -k2,2n > 5_merged_sample/${name}_OCR_fragments.bed.tmp
cd 5_merged_sample/
if [[ $# -eq 4 ]]; then
    mv ${name}_OCR_fragments.bed.tmp ${name}_OCR_fragments.bed
else
    n=`wc -l ${name}_OCR_fragments.bed.tmp | cut -f 1 -d " "`
    ratio=`bc -l <<< "${4} / $n"`
    awk -v ratio=${ratio} 'BEGIN{srand();}{if(rand() <= ratio) print $0}' ${name}_OCR_fragments.bed.tmp > ${name}_OCR_fragments.bed
    rm ${name}_OCR_fragments.bed.tmp
fi
bedToBigWig ${name}_OCR &
peak_calling ${name}_OCR ${genomeVersion}
wait
compress_bed ${name}_OCR_fragments.bed ${genomeVersion}
cd ..
for i in $@; do
   if [[ -e 2_signal/${i}_OCR_fragments.bb ]]; then
       rm 2_signal/${i}_OCR_fragments.bed
   fi
done

# ----- fragments -----
# get fragment files
IFS=',' read -r -a reads_file_array <<< ${3}
for i in ${reads_file_array[@]}; do
    if [[ ! -e 2_signal/${i}_fragments.bed ]]; then
        if [[ ! -e 2_signal/${i}_fragments.bb ]]; then
            echo "fragments file for sample ${i} do not exist, exit!"
            exit 1
        fi
        bigBedToBed 2_signal/${i}_fragments.bb 2_signal/${i}_fragments.bed
    fi
    fragments_array+=("2_signal/${i}_fragments.bed")
done

# merge fragment files
cat ${fragments_array[@]} | sort -k1,1 -k2,2n > 5_merged_sample/${name}_fragments.bed.tmp
cd 5_merged_sample/
if [[ $# -eq 4 ]]; then
    mv ${name}_fragments.bed.tmp ${name}_fragments.bed
else
    n=`wc -l ${name}_fragments.bed.tmp | cut -f 1 -d " "` && \
    ratio=`bc -l <<< "${4} / $n"` && \
    awk -v ratio=${ratio} 'BEGIN{srand();}{if(rand() <= ratio) print $0}' ${name}_fragments.bed.tmp > ${name}_fragments.bed && \
    rm ${name}_fragments.bed.tmp
fi
bedToBigWig ${name} &
peak_calling ${name} ${genomeVersion}
wait
compress_bed ${name}_fragments.bed ${genomeVersion}
cd ..
for i in $@; do
   if [[ -e 2_signal/${i}_fragments.bb ]]; then
       rm 2_signal/${i}_fragments.bed
   fi
done