#! /bin/bash

# May-28-2021

function print_help {
    echo "USAGE: $0 <name> <genomeVersion> <replicate1> <replicate2>+"
}

function bedToBigWig {
    n=`wc -l ${name}_fragments.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"`
    genomeCoverageBed -bga -scale $c -i ${name}_fragments.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_fragments.bdg
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_fragments.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}
    rm ${name}_fragments.bdg
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

MY_PATH="`dirname \"$0\"`"

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
cat ${fragments_array[@]} | sort -k1,1 -k2,2n > 4_merged_sample/${name}_fragments.bed
cd 4_merged_sample/
bedToBigWig
macs2 callpeak -f BEDPE -t ${name}_fragments.bed -n ${name} -g 1.4e9 -q 0.01 --outdir ./ --keep-dup all 2>&1 >>/dev/null | tee ${name}_MACS.out
compress_bed ${name}_fragments.bed ${genomeVersion}
cd ..

# compress fragment files
for i in $@; do
   if [[ -e 2_signal/${i}_fragments.bb ]]; then
       rm 2_signal/${i}_fragments.bed
   fi
done
