#!/usr/bin/env bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <reads1,+> <reads2,+> <main_chrom>"
    echo "Use true in 6th parameters for only map to main_chrom"
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
    name_=${1}
    genomeVersion=${2}
    genomeCoverageBed -bga -i ${name_}.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name_}.bdg && \
    ${MY_PATH}/../utilities/bdg2bw.sh ${name_}.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name_} && \
    rm ${name_}.bdg
}

# ----- parameters -----
if [[ $# -lt 6 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

name=${1}
processer=${2}
genomeVersion=${3}
reads1=${4}
reads2=${5}

case $6 in
    true ) main_chrom=true;;
    false ) main_chrom=false;;
    * ) print_help; exit 1;;
esac


MY_PATH="`readlink -f $(dirname \"$0\")`"
UTILITIES_DIR=${MY_PATH}/../utilities

IFS=',' read -r -a ReadsFiles1 <<< ${reads1}
IFS=',' read -r -a ReadsFiles2 <<< ${reads2}

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal 3_peak 4_basic_QC

trim_galore_processer=`bc <<< """${processer} / 2"""`
if [[ ${trim_galore_processer} -lt 1 ]]; then trim_galore_processer=1; fi

# ----- mapping & filtering -----
function mapping_filtering {
    echo "Step $step mapping & filtering start."
    if [[ ! -e 1_mapping/${name}.bam  ]]; then
        reads_file_process ${ReadsFiles1[@]} ${ReadsFiles2[@]}
        filteredReadsFlag=false
        for filteredReadsFile in ${filteredReads1[@]} ${filteredReads2[@]}; do
            if [[ ! -e 0_raw_data/${filteredReadsFile} ]]; then
                filteredReadsFlag=true
                break
            fi
        done
        if [[ "$filteredReadsFlag" = true ]]; then
            trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j ${trim_galore_processer} --paired ${trim_galore_input[@]} --trim-n -o 0_raw_data/ --suppress_warn
        fi
        if [[ "$main_chrom" = true ]]; then
            (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main --no-mixed  --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/Mapping_${name}.log
        else
            (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed  --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/Mapping_${name}.log
        fi
        rm ${filteredReads1[@]} ${filteredReads2[@]}
    fi
    if [[ ! -e 2_signal/${name}_raw_fragments.bed && ! -e 2_signal/${name}_raw_fragments.bb ]]; then
        bamToBed -bedpe -i 1_mapping/${name}.bam | awk '$1 !~ /_/{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}' > 2_signal/${name}_raw_fragments.bed
    fi
    if [[ ! -e 2_signal/${name}_fragments.bed && ! -e 2_signal/${name}_fragments.bb ]]; then
        sort -S 1% -k1,1 -k2,2n 2_signal/${name}_raw_fragments.bed | uniq > 2_signal/${name}_fragments.bed
    fi
    echo "Step $step end."
    step=$((step+1))
}

# ----- fragments_summary -----
function fragments_summary {
    echo "Step $step fragments summary start."
    cd 4_basic_QC
    if [[ ! -e ${name}_raw_chromosome_distribution.txt ]]; then
        ${UTILITIES_DIR}/bedCheck.sh ${name}_raw_fragments.bed && \  # assume there will be either ${name}_raw_fragments.bed or ${name}_raw_fragments.bb
        cut -f 1 ../2_signal/${name}_raw_fragments.bed | sort -S 1% | uniq -c | sort -S 1% -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${name}_raw_chromosome_distribution.txt
    fi
    if [[ ! -e ${name}_chromosome_distribution.txt ]]; then
        ${UTILITIES_DIR}/bedCheck.sh ${name}_fragments.bed  # assume there will be either ${name}_fragments.bed or ${name}_fragments.bb
        cut -f 1 ../2_signal/${name}_fragments.bed | sort -S 1% | uniq -c | sort -S 1% -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > ${name}_chromosome_distribution.txt
    fi
    if [[ ! -e ${name}_fragments_length.txt ]]; then
        ${UTILITIES_DIR}/bedCheck.sh ${name}_fragments.bed  # assume there will be either ${name}_fragments.bed or ${name}_fragments.bb
        awk '{print $3-$2}' ../2_signal/${name}_fragments.bed | sort -S 1% | uniq -c | sort -S 1% -k2,2g | awk 'BEGIN{print "fragment_length\tnumber"} {print $2"\t"$1}' > ${name}_fragments_length.txt
    fi
    cd ..
    echo "Step $step fragments summary end."
    step=$((step+1))
}

# ----- cut_sites -----
function cut_sites {
    echo "Step $step cut sites start."
    cd 2_signal
    if ([ ! -e ${name}_cutsites.bed ] && [ ! -e ${name}_cutsites.bb ]) || ([ ! -e ${name}_cutsites_plus.bed ] && [ ! -e ${name}_cutsites_plus.bb ]) || ([ ! -e ${name}_cutsites_minus.bed ] && [ ! -e ${name}_cutsites_minus.bb ]); then
        ${UTILITIES_DIR}/bedCheck.sh ${name}_fragments.bed  # assume there will be either ${name}_fragments.bed or ${name}_fragments.bb
        awk '{print $1"\t"$2"\t"$2+1"\tcut_site\t0\t+\n"$1"\t"$3-1"\t"$3"\tcut_site\t0\t-"}' ${name}_fragments.bed | sort -k1,1 -k2,2n > ${name}_cutsites.bed
        awk '{if($6=="+") print $0}' ${name}_cutsites.bed | sort -k1,1 -k2,2n > ${name}_cutsites_plus.bed
        awk '{if($6=="-") print $0}' ${name}_cutsites.bed | sort -k1,1 -k2,2n > ${name}_cutsites_minus.bed
    fi
    for name_ in ${name}_cutsites ${name}_cutsites_plus ${name}_cutsites_minus; do
        if [[ ! -e ${name_}.bw ]]; then
            bedToBigWig ${name_} ${genomeVersion}
        fi
    done
    cd ..
    echo "Step $step cut sites end."
    step=$((step+1))
}

# ----- pileup -----
function piling_up {
    echo "Step $step piling up start."
    if [[ ! -e 2_signal/${name}.bw ]]; then
        cd 2_signal
        ${UTILITIES_DIR}/bedCheck.sh ${name}_fragments.bed  # assume there will be either ${name}_fragments.bed or ${name}_fragments.bb
        fragment_length=`awk 'BEGIN{s=0;c=0} NR>1{s+=$1*$2;c+=$2} END{printf "%f", s/c}' ../4_basic_QC/${name}_fragments_length.txt`
        ${MY_PATH}/../utilities/ShiftPairEnd.sh ${name}_fragments.bed ${fragment_length} ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes && \
        n=`wc -l ${name}_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes | awk '{if($3>$2) print$0}' > ${name}_fragments_shift.bdg && \
        ${MY_PATH}/../utilities/bdg2bw.sh ${name}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_fragments_shift.bed ${name}_fragments_shift.bdg
        cd ..
    fi &
    wait
    echo "Step $step piling up end."
    step=$((step+1))
}

# ----- peak calling -----
function peak_calling {
    echo "Step $step peak calling start."
    if [[ ! -e 3_peak/${name}_MACS.out ]]; then
        cd 3_peak
        ${UTILITIES_DIR}/bedCheck.sh ${name}_fragments.bed && \  # assume there will be either ${name}_fragments.bed or ${name}_fragments.bb
        chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
        macs3 callpeak -f BEDPE -t ../2_signal/${name}_fragments.bed  --outdir ./ -n ${name} -g ${chromsize} 2>&1 >>/dev/null | tee ./${name}_MACS.out
        cd ..
    fi
    echo "Step $step peak calling end."
    step=$((step+1))
}

# ----- short fragments -----
function short_fragments {
    echo "Step $step short fragments start."
    if [[ ! -e 2_signal/${name}_OCR.bw ]]; then
        cd 2_signal
        fragment_length=`awk 'BEGIN{s=0;c=0} NR>1{if($1<=120) {s+=$1*$2;c+=$2}} END{printf "%d", s/c}' ../4_basic_QC/${name}_fragments_length.txt`
        awk '{if($3-$2<=120) print}' ${name}_fragments.bed > ${name}_OCR_fragments.bed  # assume there will be either ${name}_OCR_fragments.bed or ${name}_OCR_fragments.bb
        ${MY_PATH}/../utilities/ShiftPairEnd.sh ${name}_OCR_fragments.bed ${fragment_length} ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes
        n=`wc -l ${name}_OCR_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_OCR_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes | awk '{if($3>$2) print$0}' > ${name}_OCR_fragments_shift.bdg && \
        ${MY_PATH}/../utilities/bdg2bw.sh ${name}_OCR_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}_OCR && \
        rm ${name}_OCR_fragments_shift.bdg ${name}_OCR_fragments_shift.bed
        cd ..
    fi
    if [[ ! -e 3_peak/${name}_OCR_MACS.out ]]; then
        cd 2_signal
        ${UTILITIES_DIR}/bedCheck.sh ${name}_OCR_fragments.bed && \  # assume there will be either ${name}_OCR_fragments.bed or ${name}_OCR_fragments.bb
        chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
        macs3 callpeak -f BEDPE -t ${name}_OCR_fragments.bed  --outdir ../3_peak -n ${name}_OCR -g ${chromsize} 2>&1 >>/dev/null | tee ../3_peak/${name}_OCR_MACS.out
        cd ..
    fi
    echo "Step $step short fragments end."
    step=$((step+1))
}

# ----- clearning_up -----
function clearning_up {
    echo "Step $step clearning up start."
    rm 1_mapping/${name}.bam
    ${UTILITIES_DIR}/compressBed.sh 2_signal/${name}_raw_fragments.bed ${genomeVersion}
    ${UTILITIES_DIR}/compressBed.sh 2_signal/${name}_fragments.bed ${genomeVersion}
    ${UTILITIES_DIR}/compressBed.sh 2_signal/${name}_OCR_fragments.bed ${genomeVersion}
    echo "Step $step clearning up end."
    step=$((step+1))
}

# ----- running -----
step=1
mapping_filtering
fragments_summary
cut_sites
piling_up
peak_calling
short_fragments
clearning_up
