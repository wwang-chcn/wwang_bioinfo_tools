#! /bin/bash

# Nov-17-2020

function print_help {
        echo "USAGE: $0 <name> <genomeVersion> <replicate1> <replicate2>+"
}

function bedToBigWig {
    n=`wc -l ${name}_OCR_SE_reads.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"`
    genomeCoverageBed -bga -scale $c -i ${name}_OCR_SE_reads.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_OCR_SE_reads.bdg
    bdg2bw.sh ${name}_OCR_SE_reads.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}
    rm ${name}_OCR_SE_reads.bdg
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

if [[ $# -lt 4 ]]; then
    echo "No enought parameters!"
    print_help
    exit 1
fi

mkdir -p 3_merged_signal/

name=${1}
genomeVersion=${2}
shift 2
reads_file_array=()
for i in $@; do
    if [[ ! -e 2_signal/OCR/${i}_uniq_OCR_SE_reads.bed ]]; then
        if [[ ! -e 2_signal/OCR/${i}_uniq_OCR_SE_reads.bb ]]; then
            echo "reads for sample ${i} do not exist, exit!"
            exit 1
        fi
        bigBedToBed 2_signal/OCR/${i}_uniq_OCR_SE_reads.bb 2_signal/OCR/${i}_uniq_OCR_SE_reads.bed
    fi
    fragments_array+=("2_signal/OCR/${i}_uniq_OCR_SE_reads.bed")
done
cat ${fragments_array[@]} | sort -k1,1 -k2,2n > 3_merged_signal/${name}_OCR_SE_reads.bed
cd 3_merged_signal/
bedToBigWig
compress_bed ${name}_OCR_SE_reads.bed ${genomeVersion}
cd ..
for i in $@; do
    if [[ -e 2_signal/OCR/${i}_uniq_OCR_SE_reads.bb ]]; then
        rm 2_signal/OCR/${i}_uniq_OCR_SE_reads.bed
    fi
done