#! /bin/bash

# Dec-18-2020

function print_help {
    echo "$0 <name> <processer> <sourceDir> <genomeVersion> <reads1,+> <reads2,+>"
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

# ----- parameters -----
if [[ $# -lt 6 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi
name=${1}
processer=${2}
sourceDir=${3}
genomeVersion=${4}
reads1=${5}
reads2=${6}

MY_PATH="`dirname \"$0\"`"

IFS=',' read -r -a readsFiles1 <<< ${reads1}
IFS=',' read -r -a readsFiles2 <<< ${reads2}

mkdir -p 1_mapping 2_methylation_value 

reads_file_process ${readsFiles1[@]} ${readsFiles2[@]}
trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" --paired ${trim_galore_input[@]} --trim-n -j 4 -o 0_raw_data/ --no_report_file --suppress_warn
zcat ${filteredReads1[@]} | gzip - > 0_raw_data/${name}_1.fastq.gz &
zcat ${filteredReads2[@]} | gzip - > 0_raw_data/${name}_2.fastq.gz &
wait
bsmap -p ${processer} -a 0_raw_data/${name}_1.fastq.gz -b 0_raw_data/${name}_2.fastq.gz -d ${sourceDir}/${genomeVersion}_main.fa -R -n 1 -o 1_mapping/${name}.sam 2>&1 >>/dev/null | tee 1_mapping/Mapping_${name}.log
mcall -p ${processer} -m 1_mapping/${name}.sam --outputDir 2_methylation_value -r ${sourceDir}/${genomeVersion}_main.fa 2>&1 >>/dev/null | tee 2_methylation_value/Mcall_${name}.log && mv 1_mapping/${name}.sam.G.bed 1_mapping/${name}.sam.HG.bed 1_mapping/${name}.sam_stat.txt 2_methylation_value
rm ${filteredReads1[@]} ${filteredReads2[@]} 0_raw_data/${name}_1.fastq.gz 0_raw_data/${name}_2.fastq.gz
cd 1_mapping
samtools view -bS 1_mapping/${name}.sam > 1_mapping/${name}.bam && rm 1_mapping/${name}.sam
cd ..