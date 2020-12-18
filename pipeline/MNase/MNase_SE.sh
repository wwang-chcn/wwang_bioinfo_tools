#! /bin/bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <reads,+>"
}
function join_by {
    local IFS="$1"; shift; echo "$*"
}
function reads_file_process {
    sampleReadFiles=()
    bowtie2InputFiles=()
    for i in $@; do splitFileName=`echo ${i%.*}`; splitFileName=`echo ${splitFileName%.*}`; sampleReadFiles+=(0_raw_data/${i}); bowtie2InputFiles+=(0_raw_data/${splitFileName}_trimmed.fq.gz); done
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

IFS=',' read -r -a readsFile <<< ${reads}

mkdir -p 1_mapping 2_signal 3_basic_QC
MY_PATH="`dirname \"$0\"`"

# ----- process -----
function mapping_filtering {
    if [[ ! -e 1_mapping/${name}.bam ]]; then
        reads_file_process ${readsFile[@]} && \
        trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" --trim-n -j 4 -o 0_raw_data/ --no_report_file ${sampleReadFiles} && \
        (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed --no-unal -U ${bowtie2InputFiles} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log && \
        IFS=',' read -r -a filteredReads <<< ${bowtie2InputFiles} && \
        rm ${filteredReads[@]}
    fi
}


function piling_up {
    cd 2_signal
    if [[ ! -e ${name}.bw ]]; then
        bamToBed -i ../1_mapping/${name}.bam | awk '$1 !~ /_/{if($3>$2) print}' | sort -S 1% -k1,1 -k2,2g | uniq > ${name}_reads.bed && \
        nucleosomeShiftSE.sh ${name}_reads.bed && \
        n=`wc -l ${name}_reads_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_reads_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_reads_shift.bdg && \
        bdg2bw.sh ${name}_reads_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_reads_shift.bed ${name}_reads_shift.bdg
    fi
    cd ..
}


# ----- clearning_up -----
function clearning_up {
    rm 1_mapping/${name}.bam
    compress_bed 2_signal/${name}_reads.bed ${genomeVersion}
    #if [[ $# -eq 8 ]]; then
    #    rm 1_mapping/${name}_${SNP_strain1}.bam 1_mapping/${name}_${SNP_strain2}.bam
    #    sort -k1,1 -k2,2n 2_signal/${name}_${SNP_strain1}_fragments.bed > tmp.bed
    #    bedToBigBed -type=bed3+3 tmp.bed ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes 2_signal/${name}_${SNP_strain1}_fragments.bb
    #    rm tmp.bed 2_signal/${name}_${SNP_strain1}_fragments.bed
    #    sort -k1,1 -k2,2n 2_signal/${name}_${SNP_strain2}_fragments.bed > tmp.bed
    #    bedToBigBed -type=bed3+3 tmp.bed ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes 2_signal/${name}_${SNP_strain2}_fragments.bb
    #    rm tmp.bed 2_signal/${name}_${SNP_strain2}_fragments.bed
    #fi
}


# ----- running -----
mapping_filtering
piling_up
clearning_up
