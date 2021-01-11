#! /bin/bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <ChIPsample1> <ChIPsample2> <main_chrom> [SNP_info] [SNP_strain1] [SNP_strain2]"
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

if [[ $# -lt 6 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

name=${1}
processer=${2}
genomeVersion=${3}
ChIPsample1=${4}
ChIPsample2=${5}

case $6 in
    true ) main_chrom=true;;
    false ) main_chrom=false;;
    * ) print_help; exit 1;;
esac

IFS=',' read -r -a ChIPsampleFiles1 <<< ${ChIPsample1}
IFS=',' read -r -a ChIPsampleFiles2 <<< ${ChIPsample2}

MY_PATH="`dirname \"$0\"`"
mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal 3_peak 4_basic_QC

# ----- mapping & filtering -----
function mapping_filtering {    
    if [[ ! -e 2_signal/${name}_fragments.bed ]]; then
        if [[ ! -e 1_mapping/${name}.bam ]]; then
            reads_file_process ${ChIPsampleFiles1[@]} ${ChIPsampleFiles2[@]}
            filteredReadsFlag=false
            for filteredReadsFile in ${filteredReads1[@]} ${filteredReads2[@]}; do
                if [[ ! -e 0_raw_data/${filteredReadsFile} ]]; then
                    filteredReadsFlag=true
                    break
                fi
            done
            if [[ "$filteredReadsFlag" = true ]]; then
                trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --paired ${trim_galore_input[@]} --trim-n -o 0_raw_data/ --suppress_warn
            fi
            if [[ "$main_chrom" = true ]]; then
                (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main --no-mixed  --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log
            else
                (bowtie2 -p ${processer} --mm -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed  --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log
            fi
            rm ${filteredReads1[@]} ${filteredReads2[@]}
        fi
        bamToBed -bedpe -i 1_mapping/${name}.bam | awk '$1 !~ /_/{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3} $1 ~ /NC/{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}' > 2_signal/${name}_raw_fragments.bed
    fi
}

# ----- macs -----
function peak_calling {
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
    macs2 callpeak -f BEDPE -t 2_signal/${name}_raw_fragments.bed --outdir 3_peak -n ${name} -g ${chromsize} --nomodel --shift 37 --extsize 73 2>&1 >>/dev/null | tee 3_peak/${name}_MACS.out
}

# ----- pileup -----
function piling_up {
    if [[ ! -e 2_signal/${name}.bw ]]; then
        cd 2_signal && \
        cut -f 1-3 ${name}_raw_fragments.bed | sort -S 1% -k1,1 -k2,2n | uniq > ${name}_fragments.bed
        nucleosomeShiftPairEnd.sh ${name}_fragments.bed && \
        n=`wc -l ${name}_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name}_fragments_shift.bdg && \
        bdg2bw.sh ${name}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
        rm ${name}_fragments_shift.bed ${name}_fragments_shift.bdg
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
        bamToBed -bedpe -i ../1_mapping/${name1}.bam | awk '{if($9=="+"&&$10=="-") print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t.\t"; if($9=="-"&&$10=="+") print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t.\t"}' | awk '$1 !~ /_/{if($3>$2) print} $1 ~ /NC/{if($3>$2) print}' | sort -S 1% -k1,1 -k2,2n | uniq > ${name1}_fragments.bed && \
        nucleosomeShiftPairEnd.sh ${name1}_fragments.bed && \
        n=`wc -l ${name1}_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name1}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name1}_fragments_shift.bdg && \
        bdg2bw.sh ${name1}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name1} && \
        rm ${name1}_fragments_shift.bed ${name1}_fragments_shift.bdg
    fi &
    if [[ ! -e ${name2}.bw ]]; then
        bamToBed -bedpe -i ../1_mapping/${name2}.bam | awk '{if($9=="+"&&$10=="-") print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t.\t"; if($9=="-"&&$10=="+") print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t.\t"}' | awk '$1 !~ /_/{if($3>$2) print} $1 ~ /NC/{if($3>$2) print}' | sort -S 1% -k1,1 -k2,2n | uniq > ${name2}_fragments.bed && \
        nucleosomeShiftPairEnd.sh ${name2}_fragments.bed && \
        n=`wc -l ${name2}_fragments_shift.bed | cut -f 1 -d " "` && \
        c=`bc -l <<< "1000000 / $n"` && \
        genomeCoverageBed -bga -scale $c -i ${name2}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name2}_fragments_shift.bdg && \
        bdg2bw.sh ${name2}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name2} && \
        rm ${name2}_fragments_shift.bed ${name2}_fragments_shift.bdg
    fi &
    wait
    cd ..
}

# ----- clearning_up -----
function clearning_up {
    rm 1_mapping/${name}.bam
    compress_bed 2_signal/${name}_raw_fragments.bed ${genomeVersion}
    compress_bed 2_signal/${name}_fragments.bed ${genomeVersion}
    if [[ $# -eq 8 ]]; then
        compress_bed 2_signal/${name1}_fragments.bed ${genomeVersion}
        compress_bed 2_signal/${name2}_fragments.bed ${genomeVersion}
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