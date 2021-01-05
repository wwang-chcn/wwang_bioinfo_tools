#! /bin/bash

#! /bin/bash

function print_help {
    echo "$0 <name> <control name> <processer> <genomeVersion> <ChIPsample[,+]> <control samples[,+]>"
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

name=${1}
controlName=${2}
processer=${3}
genomeVersion=${4}
ChIPsample=${5}
ctrSample=${6}

if [[ $# -lt 6 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

IFS=',' read -r -a ChIPsampleFiles <<< ${ChIPsample}
IFS=',' read -r -a  ctrSampleFiles <<< ${ctrSample}

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal 3_peak 4_basic_QC


# ----- mapping & filtering -----
function mapping_filtering {
    if [[ ! -e 2_signal/${name}_reads.bed ]]; then
        if [[ ! -e 1_mapping/${name}.bam ]]; then
            reads_file_process ${ChIPsampleFiles[@]}
            if [[ ! -e 0_raw_data/${mapping_input_file[1]} ]]; then
                trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --trim-n -o 0_raw_data/ --suppress_warn ${trim_galore_input[@]}
            fi
            (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-unal -U ${mapping_input_file} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log && \
            rm ${filteredReads[@]} 
        fi
        bamToBed -i 1_mapping/${name}.bam | awk '$1 !~ /_/{if($3>$2) print} $1 ~ /NC/{if($3>$2) print}' > 2_signal/${name}_reads.bed
    fi
    if [[ ! -e 2_signal/${controlName}_reads.bed ]]; then
        if [[ ! -e 1_mapping/${controlName}.bam ]]; then
            reads_file_process ${ctrSampleFiles[@]}
            if [[ ! -e 0_raw_data/${mapping_input_file[1]} ]]; then
                trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --trim-n -o 0_raw_data/ --suppress_warn ${trim_galore_input[@]}
            fi
            (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-unal -U ${mapping_input_file} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${controlName}.bam) 2> 1_mapping/${controlName}_mapping.log && \
            rm ${filteredReads[@]}
        fi
        bamToBed -i 1_mapping/${controlName}.bam | awk '$1 !~ /_/{if($3>$2) print} $1 ~ /NC/{if($3>$2) print}' > 2_signal/${controlName}_reads.bed
    fi
}


# ----- macs -----
function peak_calling {
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes`
    macs2 callpeak -f BED -t 2_signal/${name}_reads.bed -c 2_signal/${controlName}_reads.bed --outdir 3_peak -n ${name} -g ${chromsize} -q 0.05 2>&1 >>/dev/null | tee 3_peak/${name}_MACS.out
}


# ----- pileup -----
function piling_up {
    cd 2_signal
    fragment_length=`grep "predicted fragment length is" ../3_peak/${name}_MACS.out | cut -f 14 -d " "`
    if [[ ! -e ${name}.bw ]]; then
        ShiftSingleEnd.sh ${name}_reads.bed ${fragment_length} && \
        n=`wc -l ${name}_reads_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_reads_shift.bdg && \
        bdg2bw.sh ${name}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_reads_shift.bed ${name}_reads_shift.bdg
    fi &
    if [[ ! -e ${controlName}.bw ]]; then
        ShiftSingleEnd.sh ${controlName}_reads.bed ${fragment_length} && \
        n=`wc -l ${controlName}_reads_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${controlName}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${controlName}_reads_shift.bdg && \
        bdg2bw.sh ${controlName}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${controlName} && \
        rm ${controlName}_reads_shift.bed ${controlName}_reads_shift.bdg
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
    compress_bed 2_signal/${name}_reads.bed ${genomeVersion}
}

# ----- running ----
mapping_filtering
peak_calling
piling_up
basic_QC
clearning_up
