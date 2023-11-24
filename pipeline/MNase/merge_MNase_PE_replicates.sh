#! /bin/bash

function print_help {
    echo "USAGE: $0 <name> <genomeVersion> <MNase replicate[,+]>"
}
function compress_bed {
    if [[ -e ${1} ]]; then
        bedFile=${1}
        col=`head -1 ${bedFile} | awk '{print NF}'`
        plus=`bc <<< "$col -3"`
        intersectBed -a ${bedFile} -b <(awk '{print $1"\t0\t"$2}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes) -wa -f 1.00 > ${bedFile}.tmp
        bedToBigBed -type=bed3+${plus} ${bedFile}.tmp ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes ${bedFile::(${#bedFile}-2)}b
        rm ${bedFile} ${bedFile}.tmp
    fi
}
function bedToBigWig {
    name_=${1}
    genomeVersion=${2}
    genomeCoverageBed -bga -i ${name_}.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name_}.bdg && \
    ${MY_PATH}/../utilities/bdg2bw.sh ${name_}.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name_} && \
    rm ${name_}.bdg
}

MY_PATH="`readlink -f $(dirname \"$0\")`"

mkdir -p 4_merged_sample/

name=${1}
genomeVersion=${2}
shift 2

# get fragment files
reads_file_array=()
for i in $@; do
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
cat ${fragments_array[@]} | sort -S 5% -k1,1 -k2,2n > 4_merged_sample/${name}_fragments.bed
cd 4_merged_sample/
bash ${MY_PATH}/../utilities/nucleosomeShiftPairEnd.sh ${name}_fragments.bed && \
n=`wc -l ${name}_fragments_shift.bed | cut -f 1 -d " "` && \
c=`bc -l <<< "1000000 / $n"` && \
genomeCoverageBed -bga -scale $c -i ${name}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_fragments_shift.bdg && \
${MY_PATH}/../utilities/bdg2bw.sh ${name}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
rm ${name}_fragments_shift.bed ${name}_fragments_shift.bdg

# cutsites from fragments file
bash ${MY_PATH}/../utilities/cutsites_from_fragments.sh ${name}_fragments.bed ${name}
for name_ in ${name}_cutsites ${name}_cutsites_plus ${name}_cutsites_minus; do
    if [[ ! -e ${name_}.bw ]]; then
        if [[ ! -e ${name_}.bed ]]; then
            bigBedToBed ${name_}.bb ${name_}.bed
        fi
        bedToBigWig ${name_} ${genomeVersion}
    fi
done

compress_bed ${name}_fragments.bed ${genomeVersion}
for name_ in ${name}_cutsites ${name}_cutsites_plus ${name}_cutsites_minus; do
    compress_bed ${name_}.bed

cd ..