#! /bin/bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <ChIPsamples,+> <main_chrom> [SNP_info] [SNP_strain1] [SNP_strain2]"
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

IFS=',' read -r -a ChIPsampleFiles <<< ${ChIPsamples}

MY_PATH="`dirname \"$0\"`"
mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal 3_peak 4_basic_QC

# ----- mapping & filtering -----
function mapping_filtering {    
    if [[ ! -e 2_signal/${name}_reads.bed ]]; then
        if [[ ! -e 1_mapping/${name}.bam ]]; then
            reads_file_process ${ChIPsampleFiles[@]}
            filteredReadsFlag=false
            for filteredReadsFile in ${filteredReads[@]}; do
                if [[ ! -e 0_raw_data/${filteredReadsFile} ]]; then
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
    macs2 callpeak -f BED -t 2_signal/${name}_raw_reads.bed --outdir 3_peak -n ${name} -g ${chromsize} --nomodel --shift 37 --extsize 73 2>&1 >>/dev/null | tee 3_peak/${name}_MACS.out
}

# ----- pileup -----
function piling_up {
    if [[ ! -e 2_signal/${name}.bw ]]; then
        cd 2_signal && \
        cut -f 1-3 ${name}_raw_reads.bed | sort -S 1% -k1,1 -k2,2n | uniq > ${name}_reads.bed
        nucleosomeShiftSE.sh ${name}_reads.bed && \
        n=`wc -l ${name}_reads_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name}_reads_shift.bdg && \
        bdg2bw.sh ${name}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_reads_shift.bed ${name}_reads_shift.bdg
        cd ..
    fi &
    wait
}

# ----- SNP -----
function SNP {
    name1=${name}_${SNP_strain1}
    name2=${name}_${SNP_strain2}
    if [[ ! -e 1_mapping/${name1}.bam || 1_mapping/${name2}.bam ]]; then
        python ${MY_PATH}/bam_split_snp.py ${SNP_info} 1_mapping/${name}.bam 1_mapping/${name1}.bam 1_mapping/${name2}.bam > 1_mapping/${name}_split_snp.log
    fi
    cd 2_signal
    if [[ ! -e ${name1}.bw ]]; then
        bamToBed -bedpe -i ../1_mapping/${name1}.bam | awk '$1 !~ /_/{print $0} $1 ~ /NC/{print $0}' | sort -S 1% -k1,1 -k2,2n | uniq > ${name1}_reads.bed && \
        nucleosomeShiftPairEnd.sh ${name1}_reads.bed && \
        n=`wc -l ${name1}_reads_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name1}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name1}_reads_shift.bdg && \
        bdg2bw.sh ${name1}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name1} && \
        rm ${name1}_reads_shift.bed ${name1}_reads_shift.bdg
    fi &
    if [[ ! -e ${name2}.bw ]]; then
        bamToBed -bedpe -i ../1_mapping/${name2}.bam| awk '$1 !~ /_/{print $0} $1 ~ /NC/{print $0}' | sort -S 1% -k1,1 -k2,2n | uniq > ${name2}_reads.bed && \
        nucleosomeShiftPairEnd.sh ${name2}_reads.bed && \
        n=`wc -l ${name2}_reads_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name2}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name2}_reads_shift.bdg && \
        bdg2bw.sh ${name2}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name2} && \
        rm ${name2}_reads_shift.bed ${name2}_reads_shift.bdg
    fi &
    wait
    cd ..
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
if [[ $# -eq 8 ]]; then
  SNP_info=${6}
  SNP_strain1=${7}
  SNP_strain2=${8}
  SNP
fi
clearning_up