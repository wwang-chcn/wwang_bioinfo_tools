#! /bin/bash

function print_help {
    echo "USAGE: $0 <name> <genomeVersion> <replicate1> <replicate2>+"
}

function bedToBigWig {
    fragment_length=`grep "predicted fragment length is" ${name}_MACS.out | cut -f 14 -d " "`
    ${MY_PATH}/../utilities/ShiftSingleEnd.sh ${name}_reads.bed ${fragment_length} ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes
    n=`wc -l ${name}_reads_shift.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"`
    genomeCoverageBed -bga -scale $c -i ${name}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_reads_shift.bdg
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}
    rm ${name}_reads_shift.bdg ${name}_reads_shift.bed
}

function compress_bed {
    bedFile=${1}
    genomeVersion=${2}
    col=`head -1 ${bedFile} | awk '{print NF}'`
    plus=`bc <<< "$col -3"`
    intersectBed -a ${bedFile} -b <(awk '{print $1"\t0\t"$2}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes) -wa -f 1.00 | sort -S 5% -k1,1 -k2,2n > ${bedFile}.tmp
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

# get read files
reads_file_array=()
for i in $@; do
    if [[ ! -e 2_signal/${i}_reads.bed ]]; then
        if [[ ! -e 2_signal/${i}_reads.bb ]]; then
            echo "fragments file for sample ${i} do not exist, exit!"
            exit 1
        fi
        bigBedToBed 2_signal/${i}_reads.bb 2_signal/${i}_reads.bed
    fi
    reads_file_array+=("2_signal/${i}_reads.bed")
done

# merge read files
cat ${reads_file_array[@]} | sort -S 5% -k1,1 -k2,2n > 4_merged_sample/${name}_reads.bed
cd 4_merged_sample/
macs2 callpeak -f BED -t ${name}_reads.bed -n ${name} -g 1.4e9 -q 0.01 --outdir ./ --keep-dup all 2>&1 >>/dev/null | tee ${name}_MACS.out
bedToBigWig
compress_bed ${name}_reads.bed ${genomeVersion}
cd ..

# compress read files
for i in $@; do
   if [[ -e 2_signal/${i}_reads.bb ]]; then
       rm 2_signal/${i}_reads.bed
   fi
done
