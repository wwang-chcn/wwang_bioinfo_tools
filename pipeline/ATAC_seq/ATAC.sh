#!/usr/bin/env bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <reads1,+> <reads2,+> [main_chrom=true]"
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
name=${1}
processer=${2}
genomeVersion=${3}
reads1=${4}
reads2=${5}

if [[ $# -lt 5 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

case $6 in
    true ) main_chrom=true;;
    false ) main_chrom=false;;
    * ) print_help; exit 1;;
esac

MY_PATH="`readlink -f $(dirname \"$0\")`"

IFS=',' read -r -a readsFiles1 <<< ${reads1}
IFS=',' read -r -a readsFiles2 <<< ${reads2}

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal

# ----- process -----
function mapping_filtering {
    if [[ ! -e 2_signal/${name}_fragments.bed ]]; then
        if [[ ! -e 1_mapping/${name}.bam ]]; then
            reads_file_process ${readsFiles1[@]} ${readsFiles2[@]}
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
                (bowtie2 -p ${processer} --trim-to 3:40 -x ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main --no-mixed --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -f 0x2 -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log
            else
                (bowtie2 -p ${processer} --trim-to 3:40 -x ~/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -f 0x2 -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/${name}_mapping.log
            fi
        rm ${filteredReads1[@]} ${filteredReads2[@]}
        fi
        bamToBed -bedpe -i 1_mapping/${name}.bam | awk '{if($9=="+"&&$10=="-") print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t."; if($9=="-"&&$10=="+") print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t."}' | awk '$1 !~ /_/{if($3>$2 && $1!="chrM") print}' | sort -S 1% -k1,1 -k2,2n > 2_signal/${name}_raw_fragments.bed &
        samtools view -@ $((${processer}-1)) -f 0x40 1_mapping/${name}.bam | cut -f 3 | sort -S 1% | uniq -c | sort -k1,1rg -S 1% | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > 2_signal/${name}_chromosome_distribution.txt &
        wait
    fi
}


function pileup { 
    cd 2_signal
    cut -f 1-3 ${name}_raw_fragments.bed | uniq > ${name}_fragments.bed
    # fragments length
    awk '{print $3-$2}' ${name}_fragments.bed | sort -S 1% | uniq -c | sort -k2,2g -S 1% | awk 'BEGIN{print "fragment_length\tnumber"} {print $2"\t"$1}' > ${name}_fragments_length.txt &
    
    # generate pesudo single-end reads
    awk 'BEGIN{srand(1007)} {if(rand()<0.5) print $1"\t"$2"\t"$2+50; else if($3-50<0) print $1"\t0\t"$3; else print $1"\t"$3-50"\t"$3}' ${name}_fragments.bed | awk '{if($2<0) print $1"\t0\t"$3; else print $0}' | sort -k1,1 -k2,2g -S 1% > ${name}_uniq_SE_reads.bed && \
    n=`wc -l ${name}_uniq_SE_reads.bed | cut -f 1 -d " "` && \
    c=`bc -l <<< "1000000 / $n"` && \
    genomeCoverageBed -bga -scale $c -i ${name}_uniq_SE_reads.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes | awk '{if($3>$2) print$0}' > ${name}_uniq_SE_reads.bdg && \
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_uniq_SE_reads.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}_uniq_SE_reads && \
    rm ${name}_uniq_SE_reads.bdg &
    wait
    # Warning: using gawk. default awk in macOS is not gawk
    cd ..
}

function OCR {
    cd 2_signal
    mkdir -p OCR
    cat ${name}_fragments.bed | awk 'BEGIN{srand(1007)} {if($3-$2<=100) {if(rand()<0.5) print $1"\t"$2"\t"$2+50"\tr"NR"\t0\t+"; else if($3-50<0) print $1"\t0\t"$3"\tr"NR"\t0\t-"; else print $1"\t"$3-50"\t"$3"\tr"NR"\t0\t-"}}' | sort -k1,1 -k2,2g > OCR/${name}_uniq_OCR_SE_reads.bed && \
    
    cd OCR && \
    n=`wc -l ${name}_uniq_OCR_SE_reads.bed | cut -f 1 -d " "` && \
    c=`bc -l <<< "1000000 / $n"` && \
    genomeCoverageBed -bga -scale $c -i ${name}_uniq_OCR_SE_reads.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes | awk '{if($3>$2) print$0}' > ${name}_uniq_OCR_SE_reads.bdg && \
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_uniq_OCR_SE_reads.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}_uniq_OCR_SE_reads && \
    rm ${name}_uniq_OCR_SE_reads.bdg &
    wait
    cd ..
}

function nucleosome {
    cd 2_signal
    mkdir -p nucleosome
    cat ${name}_fragments.bed | awk 'BEGIN{srand(1007)} {if($3-$2>=180) {if(rand()<0.5) print $1"\t"$2"\t"$2+50"\tr"NR"\t0\t+"; else if($3-50<0) print $1"\t0\t"$3"\tr"NR"\t0\t-"; else print $1"\t"$3-50"\t"$3"\tr"NR"\t0\t-"}}' | sort -k1,1 -k2,2g > nucleosome/${name}_uniq_nucleosome_SE_reads.bed && \
    
    cd nucleosome && \
    n=`wc -l ${name}_uniq_nucleosome_SE_reads.bed | cut -f 1 -d " "` && \
    c=`bc -l <<< "1000000 / $n"` && \
    genomeCoverageBed -bga -scale $c -i ${name}_uniq_nucleosome_SE_reads.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes | awk '{if($3>$2) print$0}' > ${name}_uniq_nucleosome_SE_reads.bdg && \
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_uniq_nucleosome_SE_reads.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}_uniq_nucleosome_SE_reads && \
    rm ${name}_uniq_nucleosome_SE_reads.bdg &
    wait
    cd ..
}


# ----- clearning_up -----
function clearning_up {
    rm 1_mapping/${name}.bam
    compress_bed 2_signal/${name}_fragments.bed ${genomeVersion}
    compress_bed 2_signal/${name}_raw_fragments.bed ${genomeVersion}
    compress_bed 2_signal/${name}_uniq_SE_reads.bed ${genomeVersion}
    compress_bed 2_signal/OCR/${name}_uniq_OCR_SE_reads.bed ${genomeVersion}
    compress_bed 2_signal/nucleosome/${name}_uniq_nucleosome_SE_reads.bed ${genomeVersion}
}


# ----- running ----
mapping_filtering
pileup
OCR
nucleosome
clearning_up