#! /bin/bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <ChIPsamples,+> <main_chrom>"
    echo "Use true in 5th parameters for only map to main_chrom"
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

if [[ $# -lt 5 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

name=${1}
processer=${2}
genomeVersion=${3}
ChIPsamples=${4}

case $5 in
    true ) main_chrom=true;;
    false ) main_chrom=false;;
    * ) print_help; exit 1;;
esac

MY_PATH="`readlink -f $(dirname \"$0\")`"

IFS=',' read -r -a ChIPsampleFiles <<< ${ChIPsamples}

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal 3_peak 4_basic_QC

# ----- mapping & filtering -----
function mapping_filtering {    
    if [[ ! -e 2_signal/${name}_raw_reads.bed && ! -e 2_signal/${name}_raw_reads.bb ]]; then
        if [[ ! -e 1_mapping/${name}.bam ]]; then
            reads_file_process ${ChIPsampleFiles[@]}
            filteredReadsFlag=false
            for filteredReadsFile in ${filteredReads[@]}; do
                if [[ ! -e ${filteredReadsFile} ]]; then
                    filteredReadsFlag=true
                    break
                fi
            done
            if [[ "$filteredReadsFlag" = true ]]; then
                trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --trim-n -o 0_raw_data/ --suppress_warn ${trim_galore_input[@]} 2> 0_raw_data/Trim_galore_${name}.log
            fi
            if [[ "$main_chrom" = true ]]; then
                (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main --no-mixed  --no-discordant --no-unal -U ${mapping_input_file} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log
            else
                (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed  --no-discordant --no-unal -U ${mapping_input_file} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log
            fi
            rm ${filteredReads[@]}
        fi
        bamToBed -i 1_mapping/${name}.bam | awk '$1 !~ /_/{print $0} $1 ~ /NC/{print $0}' > 2_signal/${name}_raw_reads.bed
    fi
}

# ----- macs -----
function peak_calling {
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
    macs3 callpeak -f BED -t 2_signal/${name}_raw_reads.bed --outdir 3_peak -n ${name} -g ${chromsize} --nomodel --shift 37 --extsize 73 2>&1 >>/dev/null | tee 3_peak/${name}_MACS.out
}

# ----- pileup -----
function piling_up {
    if [[ ! -e 2_signal/${name}.bw ]]; then
        cd 2_signal && \
        awk '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$6}' ${name}_raw_reads.bed | sort -S 1% -k1,1 -k2,2n | uniq | awk '{print $1"\t"$2"\t"$3"\tReads"NR"\t"$5"\t"$6}' > ${name}_reads.bed
        ${MY_PATH}/../utilities/nucleosomeShiftSingleEnd.sh ${name}_reads.bed && \
        n=`wc -l ${name}_reads_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name}_reads_shift.bdg && \
        ${MY_PATH}/../utilities/bdg2bw.sh ${name}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_reads_shift.bed ${name}_reads_shift.bdg
        cd ..
    fi &
    wait
}

# ----- clearning_up -----
function clearning_up {
    rm 1_mapping/${name}.bam
    compress_bed 2_signal/${name}_raw_reads.bed ${genomeVersion}
    compress_bed 2_signal/${name}_reads.bed ${genomeVersion}
    if [[ $# -eq 8 ]]; then
        compress_bed 2_signal/${name1}_reads.bed ${genomeVersion}
        compress_bed 2_signal/${name2}_reads.bed ${genomeVersion}
    fi
}

# ----- running -----
mapping_filtering
peak_calling
piling_up
clearning_up
