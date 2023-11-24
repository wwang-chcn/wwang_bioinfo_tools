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

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal

MY_PATH="`readlink -f $(dirname \"$0\")`"

IFS=',' read -r -a readsFiles <<< ${reads}

# ----- process -----
function mapping {
    if [[ ! -e 1_mapping/${name}.bam ]]; then
        reads_file_process ${readsFiles[@]} && \
        trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --trim-n -o 0_raw_data/ --no_report_file ${trim_galore_input[@]} 2> 0_raw_data/Trim_galore_${name}.log && \
        hisat2 -5 1 --summary-file 1_mapping/Mapping_${name}.log --new-summary -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main -U ${mapping_input_file} | samtools view -@ $((${processer}-1)) -bS > 1_mapping/${name}.bam && \
        rm ${filteredReads[@]}
    fi
}


function signal {
    if [[ ! -e 2_signal/${name}.bw ]]; then
        cd 2_signal && \
        bamToBed -i ../1_mapping/${name}.bam | awk '{if($6=="+") print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6; else print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}' > ${name}.bed && \
        sort -S 1% -k1,1 -k2,2n ${name}.bed | genomeCoverageBed -bga -i - -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}.bdg && \
        ${MY_PATH}/../utilities/bdg2bw.sh ${name}.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}.bdg
        cd ..
    fi
}



function cleaning_up {
    rm 1_mapping/${name}.bam
    compress_bed 2_signal/${name}.bed ${genomeVersion}
}


mapping
signal
cleaning_up
