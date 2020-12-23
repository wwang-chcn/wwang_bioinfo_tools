#! /bin/bash

function print_help {
    echo "$0 <name> <control name> <processer> <genomeVersion> <ChIPsample1[,+]> <ChIPsample2[,+]> <ctrsample1[,+]> <ctrsample2[,+]>"
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

name=${1}
controlName=${2}
processer=${3}
genomeVersion=${4}
ChIPsample1=${5}
ChIPsample2=${6}
ctrsample1=${7}
ctrsample2=${8}

if [[ $# -lt 8 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

IFS=',' read -r -a ChIPsampleFiles1 <<< ${ChIPsample1}
IFS=',' read -r -a ChIPsampleFiles2 <<< ${ChIPsample2}
IFS=',' read -r -a  ctrsampleFiles1 <<< ${ctrsample1}
IFS=',' read -r -a  ctrsampleFiles2 <<< ${ctrsample2}


mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal 3_peak 4_basic_QC


# ----- mapping & filtering -----
function mapping_filtering {
    if [[ ! -e 2_signal/${controlName}_fragments.bed ]]; then
        reads_file_process ${ctrsampleFiles1[@]} ${ctrsampleFiles2[@]} && \
        trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --paired ${trim_galore_input[@]} --trim-n -o 0_raw_data/ --no_report_file --suppress_warn && \
        (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${controlName}.bam) 2> 1_mapping/${controlName}_mapping.log && \
        rm ${filteredReads1[@]} ${filteredReads2[@]} && \
        bamToBed -bedpe -i 1_mapping/${controlName}.bam | awk '$1 !~ /_/{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}'  | sort -S 1% -k1,1 -k2,2n | uniq > 2_signal/${controlName}_fragments.bed
    fi
    
    if [[ ! -e 2_signal/${name}_fragments.bed ]]; then
        reads_file_process ${ChIPsampleFiles1[@]} ${ChIPsampleFiles2[@]} && \
        trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --paired ${trim_galore_input[@]} --trim-n -o 0_raw_data/ --no_report_file --suppress_warn && \
        (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed  --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log && \
        rm ${filteredReads1[@]} ${filteredReads2[@]} && \
        bamToBed -bedpe -i 1_mapping/${name}.bam | awk '$1 !~ /_/{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}'  | sort -S 1% -k1,1 -k2,2n | uniq > 2_signal/${name}_fragments.bed
    fi
}


# ----- macs -----
function peak_calling {
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes`
    macs2 callpeak -f BEDPE -t 2_signal/${name}_fragments.bed -c 2_signal/${controlName}_fragments.bed --outdir 3_peak -n ${name} -g ${chromsize} -q 0.05 2>&1 >>/dev/null | tee 3_peak/${name}_MACS.log
}


# ----- pileup -----
function piling_up {
    cd 2_signal
    fragment_length=`grep "fragment size =" ../3_peak/${name}_MACS.log | cut -f 12 -d " "`
    if [[ ! -e ${controlName}.bw ]]; then
        ShiftPairEnd.sh ${controlName}_fragments.bed ${fragment_length} && \
        n=`wc -l ${controlName}_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${controlName}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${controlName}_fragments_shift.bdg && \
        bdg2bw.sh ${controlName}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${controlName} && \
        rm ${controlName}_fragments_shift.bed ${controlName}_fragments_shift.bdg
    fi &

    if [[ ! -e ${name}.bw ]]; then
        ShiftPairEnd.sh ${name}_fragments.bed ${fragment_length} && \
        n=`wc -l ${name}_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_fragments_shift.bdg && \
        bdg2bw.sh ${name}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_fragments_shift.bed ${name}_fragments_shift.bdg
    fi &
    wait
    cd ..
}

# ----- basic_QC -----
function basic_QC {
    cd 4_basic_QC
    cut -f 1-3 ../3_peak/${name}_peaks.narrowPeak > ../3_peak/${name}_peaks.bed
    ceasBW -b ../3_peak/${name}_peaks.bed -w ../2_signal/${name}.bw -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.refGene.sq3 --name ${name} -l ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes --bg
    cd ..
}

# ----- clearning_up -----
function clearning_up {
    rm 1_mapping/${name}.bam
    compress_bed 2_signal/${name}_fragments.bed ${genomeVersion}
}


# ----- running ----
mapping_filtering
peak_calling
piling_up
basic_QC
clearning_up