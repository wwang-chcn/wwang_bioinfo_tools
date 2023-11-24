#! /bin/bash

# May-28-2021

function print_help {
    echo "USAGE: $0 <name> <genomeVersion> <replicate1> <replicate2>+"
}

function bedToBigWig {
    n=`wc -l ${name}_reads.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"`
    genomeCoverageBed -bga -scale $c -i ${name}_reads.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_reads.bdg
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_reads.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}
    rm ${name}_reads.bdg
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

if [[ $# -lt 3 ]]; then
    echo "No enought parameters!"
    print_help
    exit 1
fi

MY_PATH="`readlink -f $(dirname \"$0\")`"

mkdir -p 5_merged_sample/

name=${1}
genomeVersion=${2}
shift 2

# get reads files
reads_file_array=()
for i in $@; do
    if [[ ! -e 2_signal/${i}_reads.bed ]]; then
        if [[ ! -e 2_signal/${i}_reads.bb ]]; then
            echo "fragments file for sample ${i} do not exist, exit!"
            exit 1
        fi
        bigBedToBed 2_signal/${i}_reads.bb 2_signal/${i}_reads.bed
    fi
    fragments_array+=("2_signal/${i}_reads.bed")
done

# merge reads files
cat ${fragments_array[@]} | sort -k1,1 -k2,2n > 5_merged_sample/${name}_reads.bed
cd 5_merged_sample/
bedToBigWig
compress_bed ${name}_reads.bed ${genomeVersion}
cd ..
#for i in $@; do
#    if [[ -e 2_signal/${i}_reads.bb ]]; then
#        rm 2_signal/${i}_reads.bed
#    fi
#done
