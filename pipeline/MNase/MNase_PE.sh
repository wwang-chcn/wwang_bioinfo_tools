#! /bin/bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <reads1,+> <reads2,+> [SNP_info] [SNP_strain1] [SNP_strain2]"
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


IFS=',' read -r -a readsFiles1 <<< ${reads1}
IFS=',' read -r -a readsFiles2 <<< ${reads2}

mkdir -p 1_mapping 2_signal 3_basic_QC
MY_PATH="`dirname \"$0\"`"
SNP_info=`readlink -f ${SNP_info}`

# ----- process -----
function mapping_filtering {
    if [[ ! -e 1_mapping/${name}.bam ]]; then
        reads_file_process ${readsFiles1[@]} ${readsFiles2[@]} && \
        trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" --paired ${trim_galore_input[@]} --trim-n -j 4 -o 0_raw_data/ --no_report_file --suppress_warn && \
        (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> ${name}_mapping.log && \
        rm ${filteredReads1[@]} ${filteredReads2[@]}
    fi
}


function piling_up {
    cd 2_signal
    if [[ ! -e ${name}.bw ]]; then
        bamToBed -bedpe -i ../1_mapping/${name}.bam | awk '{if($9=="+"&&$10=="-") print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t.\t"; if($9=="-"&&$10=="+") print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t.\t"}' | awk '$1 !~ /_/{if($3>$2) print}' | sort -S 1% -k1,1 -k2,2g | uniq > ${name}_fragments.bed && \
        nucleosomeShiftPairEnd.sh ${name}_fragments.bed && \
        n=`wc -l ${name}_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_fragments_shift.bdg && \
        bdg2bw.sh ${name}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_fragments_shift.bed ${name}_fragments_shift.bdg
    fi
    cd ..
}


function SNP {
    cd 1_mapping && \
    python ../${MY_PATH}/bam_split_snp.py ${SNP_info} ${name}.bam ${name}_${SNP_strain1}.bam ${name}_${SNP_strain2}.bam > ${name}_split_snp.log && \
    cd ../2_signal && \
        if [[ ! -e ${name}_${SNP_strain1}.bw ]]; then
            bamToBed -bedpe -i ../1_mapping/${name}_${SNP_strain1}.bam | awk '{if($9=="+"&&$10=="-") print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t.\t"; if($9=="-"&&$10=="+") print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t.\t"}' | awk '$1 !~ /_/{if($3>$2) print}' | sort -S 1% -k1,1 -k2,2g | uniq > ${name}_${SNP_strain1}_fragments.bed && \
            nucleosomeShiftPairEnd.sh ${name}_${SNP_strain1}_fragments.bed && \
            n=`wc -l ${name}_${SNP_strain1}_fragments_shift.bed | cut -f 1 -d " "` && \
            c=`bc -l <<< "1000000 / $n"` && \
            genomeCoverageBed -bga -scale $c -i ${name}_${SNP_strain1}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_${SNP_strain1}_fragments_shift.bdg && \
            bdg2bw.sh ${name}_${SNP_strain1}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}_${SNP_strain1} && \
            rm ${name}_${SNP_strain1}_fragments_shift.bed ${name}_${SNP_strain1}_fragments_shift.bdg
        fi
        if [[ ! -e ${name}_${SNP_strain2}.bw ]]; then
            bamToBed -bedpe -i ../1_mapping/${name}_${SNP_strain2}.bam | awk '{if($9=="+"&&$10=="-") print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t.\t"; if($9=="-"&&$10=="+") print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t.\t"}' | awk '$1 !~ /_/{if($3>$2) print}' | sort -S 1% -k1,1 -k2,2g | uniq > ${name}_${SNP_strain2}fragments.bed && \
            nucleosomeShiftPairEnd.sh ${name}_${SNP_strain2}_fragments.bed && \
            n=`wc -l ${name}_${SNP_strain2}_fragments_shift.bed | cut -f 1 -d " "` && \
            c=`bc -l <<< "1000000 / $n"` && \
            genomeCoverageBed -bga -scale $c -i ${name}_${SNP_strain2}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_${SNP_strain2}_fragments_shift.bdg && \
            bdg2bw.sh ${name}_${SNP_strain2}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}_${SNP_strain2} && \
            rm ${name}_${SNP_strain2}_fragments_shift.bed ${name}_${SNP_strain2}_fragments_shift.bdg
        fi
    cd ..
}

# ----- clearning_up -----
function clearning_up {
    rm 1_mapping/${name}.bam
    compress_bed 2_signal/${name}_fragments.bed ${genomeVersion}
    if [[ $# -eq 8 ]]; then
        rm 1_mapping/${name}_${SNP_strain1}.bam 1_mapping/${name}_${SNP_strain2}.bam
        compress_bed 2_signal/${name}_${SNP_strain1}_fragments.bed ${genomeVersion}
        compress_bed 2_signal/${name}_${SNP_strain2}_fragments.bed ${genomeVersion}
    fi
}

# ----- running -----
mapping_filtering
piling_up
if [[ $# -eq 8 ]]; then
    SNP_info=${6}
    SNP_strain1=${7}
    SNP_strain2=${8}
    SNP
fi
clearning_up

