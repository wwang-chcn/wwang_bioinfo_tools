#! /bin/bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <R1.fastq.gz,+> <R2.fastq.gz,+> [SNP_info] [strain1] [strain2]"
}
function join_by {
    local IFS="$1"; shift; echo "$*"
}
function reads_file_process {
    sampleReadFiles1=()
    sampleReadFiles2=()
    for i in `seq $(($#/2))`; do sampleReadFiles1+=(${@:$i:1}); done
    for i in `seq $(($#/2+1)) $#`; do sampleReadFiles2+=(${@:$i:1}); done
    trim_galore_input=()
    for index in ${!sampleReadFiles1[@]}; do trim_galore_input+=(0_raw_data/${sampleReadFiles1[$index]} 0_raw_data/${sampleReadFiles2[$index]}); done
    filteredReads1=()
    for readFile in ${sampleReadFiles1[@]}; do splitFileName=`echo ${readFile%.*}`; splitFileName=`echo ${splitFileName%.*}`; filteredReads1+=(0_raw_data/${splitFileName}_val_1.fq.gz); done
    mapping_input_file1=`join_by , ${filteredReads1[@]}`
    filteredReads2=()
    for readFile in ${sampleReadFiles2[@]}; do splitFileName=`echo ${readFile%.*}`; splitFileName=`echo ${splitFileName%.*}`; filteredReads2+=(0_raw_data/${splitFileName}_val_2.fq.gz); done
    mapping_input_file2=`join_by , ${filteredReads2[@]}`
}
function bedToBigWig {
    n=`cut -f 4 ${3}.bed | sort -S 1% | uniq | wc -l` && \
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
if [[ $# -lt 5 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

name=${1}
processer=${2}
genomeVersion=${3}
reads1=${4}
reads2=${5}

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_expression_value 3_signal 4_basic_QC
MY_PATH="`dirname \"$0\"`"

IFS=',' read -r -a readsFiles1 <<< ${reads1}
IFS=',' read -r -a readsFiles2 <<< ${reads2}

# ----- process -----
function mapping {
    if [[ ! -e 1_mapping/${name}.bam ]]; then
        reads_file_process ${readsFiles1[@]} ${readsFiles2[@]} && \
        trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" --paired ${trim_galore_input[@]} --trim-n -j 4 -o 0_raw_data/ --no_report_file --suppress_warn && \
        hisat2 --summary-file 1_mapping/Mapping_${name}.log --new-summary -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bS > 1_mapping/${name}.bam && \
        rm ${filteredReads1[@]} ${filteredReads2[@]}
    fi
}

function expression_value {
    cd 1_mapping/ && \
    samtools sort -@ $((${processer}-1)) -o ${name}_sorted.bam ${name}.bam && samtools index -@ $((${processer}-1)) ${name}_sorted.bam && \
    cd ../2_expression_value
    stringtie ../1_mapping/${name}_sorted.bam -p $((${processer}/2)) -G ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.gtf -l ${name} -A RNA_seq_${name}_refGene_coverage.txt -o RNA_seq_${name}_refGene_coverage.gtf -e -B &
    stringtie ../1_mapping/${name}_sorted.bam -p $((${processer}/2)) -G ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.gtf -l ${name} -A RNA_seq_${name}_ensGene_coverage.txt -o RNA_seq_${name}_ensGene_coverage.gtf -e -B &
    samtools view ../1_mapping/${name}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.gtf  -tag stdin -o ${name}_refGene.read_cnt &
    samtools view ../1_mapping/${name}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.gtf  -tag stdin -o ${name}_ensGene.read_cnt &
    wait
    rm ../1_mapping/${name}_sorted.bam ../1_mapping/${name}_sorted.bam.bai
    cd ..
}

function signal {
    if [[ ! -e 3_signal/RNA_seq_${name}.bw ]]; then
        cd 3_signal && \
        cat ../1_mapping/${name}.bam | python ~/bin/myscripts/bamToBed.py > RNA_seq_${name}.bed && \
        bedToBigWig RNA_seq_${name} RNA_seq_${name} RNA_seq_${name} && \
        cd ..
    fi
}

function basic_QC {
    cd 4_basic_QC
    read_distribution.py -i ../1_mapping/${name}.bam -r ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.bed > RNA_seq_${name}_refGene_read_distribution.txt &
    read_distribution.py -i ../1_mapping/${name}.bam -r ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.bed > RNA_seq_${name}_ensGene_read_distribution.txt &
    geneBody_coverage2.py -i ../3_signal/RNA_seq_${name}.bw -r ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.bed -o RNA_seq_${name}_refGene &
    geneBody_coverage2.py -i ../3_signal/RNA_seq_${name}.bw -r ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.bed -o RNA_seq_${name}_ensGene &
    wait
    cd ..
}

function SNP {
   name1=${name}_${SNP_strain1}
   name2=${name}_${SNP_strain2}
   if [[ ! -e 1_mapping/${name1}.bam || ! -e 1_mapping/${name2}.bam ]]; then
        python ${MY_PATH}/bam_split_snp.py ${SNP_info} 1_mapping/${name}.bam 1_mapping/${name1}.bam 1_mapping/${name2}.bam > 1_mapping/${name}_split_snp.log
   fi
   if [[ ! -e 3_signal/RNA_seq_${name1}.bw ]]; then
        cd 1_mapping/ && \
        samtools sort -@ $((${processer}-1)) -o ${name1}_sorted.bam ${name1}.bam && samtools index -@ $((${processer}-1)) ${name1}_sorted.bam && \
        cd ../2_expression_value && \
        stringtie ../1_mapping/${name1}_sorted.bam -p $((${processer}/2)) -G ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.gtf -l ${name1} -A ${name1}_refGene_coverage.txt -o ${name1}_refGene_coverage.gtf -e -B && \
        stringtie ../1_mapping/${name1}_sorted.bam -p $((${processer}/2)) -G ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.gtf -l ${name1} -A ${name1}_ensGene_coverage.txt -o ${name1}_ensGene_coverage.gtf -e -B && \
        rm ../1_mapping/${name1}_sorted.bam ../1_mapping/${name1}_sorted.bam.bai && \
        samtools view ../1_mapping/${name1}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.gtf  -tag stdin -o ${name1}_refGene.read_cnt && \
        samtools view ../1_mapping/${name1}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.gtf  -tag stdin -o ${name1}_ensGene.read_cnt && \
        cd ../3_signal && \
        cat ../1_mapping/${name1}.bam | python ~/bin/myscripts/bamToBed.py | awk '$1 !~ /_/{print}' > RNA_seq_${name1}.bed && \
        bedToBigWig RNA_seq_${name1} RNA_seq_${name1} RNA_seq_${name1}
        cd ..
    fi &
    if [[ ! -e 3_signal/RNA_seq_${name2}.bw ]]; then
        cd 1_mapping/ && \
        samtools sort -@ $((${processer}-1)) -o ${name2}_sorted.bam ${name2}.bam && samtools index -@ $((${processer}-1)) ${name2}_sorted.bam && \
        cd ../2_expression_value && \
        stringtie ../1_mapping/${name2}_sorted.bam -p $((${processer}/2)) -G ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.gtf -l ${name2} -A ${name2}_refGene_coverage.txt -o ${name2}_refGene_coverage.gtf -e -B && \
        stringtie ../1_mapping/${name2}_sorted.bam -p $((${processer}/2)) -G ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.gtf -l ${name2} -A ${name2}_ensGene_coverage.txt -o ${name2}_ensGene_coverage.gtf -e -B && \
        rm ../1_mapping/${name2}_sorted.bam ../1_mapping/${name2}_sorted.bam.bai && \
        samtools view ../1_mapping/${name2}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.gtf  -tag stdin -o ${name2}_refGene.read_cnt && \
        samtools view ../1_mapping/${name2}.bam | gfold count -ann ~/source/bySpecies/${genomeVersion}/${genomeVersion}.ensGene.gtf  -tag stdin -o ${name2}_ensGene.read_cnt && \
        cd ../3_signal && \
        cat ../1_mapping/${name2}.bam | python ~/bin/myscripts/bamToBed.py | awk '$1 !~ /_/{print}' > RNA_seq_${name2}.bed && \
        bedToBigWig RNA_seq_${name2}_fragments RNA_seq_${name2} RNA_seq_${name2}
        cd ..
    fi &
    wait
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
if [[ $# -eq 8 ]]; then
    SNP_info=${6}
    SNP_strain1=${7}
    SNP_strain2=${8}
    SNP
fi
cleaning_up