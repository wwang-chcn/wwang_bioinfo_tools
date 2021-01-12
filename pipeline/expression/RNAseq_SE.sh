#! /bin/bash

# May-30-2019

# ----- functions -----
function print_help {
    echo "$0 <name> <processer> <genomeVersion> <readsFile+>"
}

function join_by {
    local IFS="$1"; shift; echo "$*"
}

function reads_file_process {
    sampleReadFiles=($@)
    trim_galore_input=()
    for index in ${!sampleReadFiles[@]}; do trim_galore_input+=(0_raw_data/${sampleReadFiles[$index]}); done
    filteredReads=()
    for readFile in ${sampleReadFiles[@]}; do splitFileName=`echo ${readFile%.*}`; splitFileName=`echo ${splitFileName%.*}`; filteredReads+=(0_raw_data/${splitFileName}_trimmed.fq.gz); done
    mapping_input_file=`join_by , ${filteredReads[@]}`
}

function bedToBigWig {
    n=`cut -f 4 ${1}.bed | sort -S 1% | uniq | wc -l` && \
    c=`bc -l <<< "1000000 / ${n}"` && \
    sort -S 1% -k1,1 -k2,2n ${1}.bed | genomeCoverageBed -bga -scale $c -i - -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${1}.bdg && \
    /mnt/Storage/home/wangwen/bin/myscripts/bdg2bw.sh ${1}.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${2} && \
    rm ${1}.bdg
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

# ----- parameters -----
if [[ $# -lt 4 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

name=${1}
processer=${2}
genomeVersion=${3}
reads=${4}

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_expression_value 3_signal 4_basic_QC
MY_PATH="`dirname \"$0\"`"

IFS=',' read -r -a readsFiles <<< ${reads}

# ----- process -----
function mapping {
    if [[ ! -e 1_mapping/${name}.bam ]]; then
        reads_file_process ${readsFiles[@]} && \
        trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --trim-n -o 0_raw_data/ --no_report_file ${trim_galore_input[@]} 2> 0_raw_data/Trim_galore_${name}.log && \
        hisat2 --summary-file 1_mapping/Mapping_${name}.log --new-summary -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main -U ${mapping_input_file} | samtools view -@ $((${processer}-1)) -bS > 1_mapping/${name}.bam && \
        rm ${filteredReads[@]}
    fi
}

function expression_value {
    cd 1_mapping/ && \
    samtools sort -@ $((${processer}-1)) -o ${name}_sorted.bam ${name}.bam && samtools index ${name}_sorted.bam && \
    cd ../2_expression_value
    stringtie ../1_mapping/${name}_sorted.bam -p $((${processer}/2)) -G ~/source/Annotation/zv9.ens.gtf -l ${name} -A RNA_seq_${name}_ensGene_coverage.txt -o RNA_seq_${name}_ensGene_coverage.gtf -e -B &
    stringtie ../1_mapping/${name}_sorted.bam -p $((${processer}/2)) -G ~/source/Annotation/RNASeq_56535_embryonic_transcriptome_merged.gtf -l ${name} -A RNA_seq_${name}_RNAanno_coverage.txt -o RNA_seq_${name}_RNAanno_coverage.gtf -e -B &
    samtools view ../1_mapping/${name}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.gtf  -tag stdin -o ${name}_ensGene.read_cnt &
    samtools view ../1_mapping/${name}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/RNASeq_56535_embryonic_transcriptome_merged.gtf  -tag stdin -o ${name}_RNAanno.read_cnt &
    wait
    rm ../1_mapping/${name}_sorted.bam ../1_mapping/${name}_sorted.bam.bai
    cd ..
}

function signal {
    if [[ ! -e 3_signal/RNA_seq_${name}.bw ]]; then
        cd 3_signal && \
        cat ../1_mapping/${name}.bam | python ~/bin/myscripts/bamToBed.py > RNA_seq_${name}.bed && \
        bedToBigWig RNA_seq_${name} RNA_seq_${name} && \
        cd ..
    fi
}

function basic_QC {
    cd 4_basic_QC
    read_distribution.py -i ../1_mapping/${name}.bam -r                       ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.bed > RNA_seq_${name}_ensGene_read_distribution.txt &
    read_distribution.py -i ../1_mapping/${name}.bam -r ~/source/bySpecies/${genomeVersion}/RNASeq_merged_embryonic_transcriptome.bed > RNA_seq_${name}_RNAanno_read_distribution.txt &
    geneBody_coverage2.py -i ../3_signal/RNA_seq_${name}.bw -r                       ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.bed -o RNA_seq_${name}_ensGene &
    geneBody_coverage2.py -i ../3_signal/RNA_seq_${name}.bw -r ~/source/bySpecies/${genomeVersion}/RNASeq_merged_embryonic_transcriptome.bed -o RNA_seq_${name}_RNAanno &
    #wait
    cd ..
}

function cleaning_up {
    rm 1_mapping/${name}.bam
    compress_bed 3_signal/RNA_seq_${name}.bed ${genomeVersion}
    if [[ $# -eq 8 ]]; then
        rm 1_mapping/${name1}.bam 1_mapping/${name2}.bam
        compress_bed 3_signal/RNA_seq_${name1}.bed ${genomeVersion}
        compress_bed 3_signal/RNA_seq_${name2}.bed ${genomeVersion}
    fi
}


mapping
expression_value
signal
basic_QC
cleaning_up
