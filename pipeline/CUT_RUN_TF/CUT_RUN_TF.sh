#!/usr/bin/bash

function print_help {
    echo "$0 <name> <processer> <genomeVersion> <reads1,+> <reads2,+> <main_chrom> [] [] []"
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

IFS=',' read -r -a ReadsFiles1 <<< ${reads1}
IFS=',' read -r -a ReadsFiles2 <<< ${reads2}

if [[ $# -eq 9 ]]; then
    SNP_info=${7}
    SNP_strain1=${8}
    SNP_strain2=${9}
    name1=${name}_${SNP_strain1}
    name2=${name}_${SNP_strain2}
fi

mkdir -p 0_raw_data/FastQC_OUT 1_mapping 2_signal 3_peak
MY_PATH="`dirname \"$0\"`"


# ----- mapping & filtering -----

function mapping_filtering {
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
            trim_galore --fastqc --fastqc_args "--outdir 0_raw_data/FastQC_OUT --nogroup -t ${processer} -q" -j 4 --paired ${trim_galore_input[@]} --trim-n -o 0_raw_data/ --suppress_warn
        fi
        if [[ "$main_chrom" = true ]]; then
            (bowtie2 -p ${processer} --mm -x /mnt/Storage/home/wangwen/source/bySpecies/${genomeVersion}/${genomeVersion}_main --no-mixed  --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/Mapping_${name}.log
        else
            (bowtie2 -p ${processer} --mm -x /mnt/Storage/home/wangwen/source/bySpecies/${genomeVersion}/${genomeVersion} --no-mixed  --no-discordant --no-unal -1 ${mapping_input_file1} -2 ${mapping_input_file2} | samtools view -@ $((${processer}-1)) -bSq 30 > 1_mapping/${name}.bam) 2> 1_mapping/Mapping_${name}.log
        fi
        rm ${filteredReads1[@]} ${filteredReads2[@]}
    fi
    if [[ ! -e 2_signal/${name}_fragments.bed  ]]; then
      bamToBed -bedpe -i 1_mapping/${name}.bam | awk '$1 !~ /_/{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}'  | sort -S 1% -k1,1 -k2,2n | uniq > 2_signal/${name}_fragments.bed
      cut -f 1 2_signal/${name}_fragments.bed | sort -S 1% | uniq -c | sort -S 1% -k1,1rg | awk 'BEGIN{print "chromosome\tnumber"} {print $2"\t"$1}' > 2_signal/${name}_chromosome_distribution.txt
      awk '{print $3-$2}' 2_signal/${name}_fragments.bed | sort -S 1% | uniq -c | sort -S 1% -k2,2g | awk 'BEGIN{print "fragment_length\tnumber"} {print $2"\t"$1}' > 2_signal/${name}_fragments_length.txt
    fi
}


# ----- pileup -----
function piling_up {
  cd 2_signal
  fragment_length=`awk 'BEGIN{s=0;c=0} NR>1{s+=$1*$2;c+=$2} END{printf "%f", s/c}' ${name}_fragments_length.txt`
  if [[ ! -e ${name}.bw ]]; then
    ShiftPairEnd.sh ${name}_fragments.bed ${fragment_length} && \
    n=`wc -l ${name}_fragments_shift.bed | cut -f 1 -d " "` && \
    c=`bc -l <<< "1000000 / $n"` && \
    genomeCoverageBed -bga -scale $c -i ${name}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes > ${name}_fragments_shift.bdg && \
    bdg2bw.sh ${name}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
    rm ${name}_fragments_shift.bed ${name}_fragments_shift.bdg
  fi &
  wait
  cd ..
}


# ----- macs -----
function peak_calling {
  chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes`
  macs2 callpeak -f BEDPE -t 2_signal/${name}_fragments.bed --outdir 3_peak -n ${name} -g ${chromsize} -q 0.05 2>&1 >>/dev/null | tee 2_peak/MACS_${name}.log
}


# ----- short fragments -----
function short_fragments {

  cd 2_signal
  fragment_length=`awk 'BEGIN{s=0;c=0} NR>1{if($1<=120) {s+=$1*$2;c+=$2}} END{printf "%d", s/c}' ${name}_fragments_length.txt`
  awk '{if($3-$2<=120) print}' ${name}_fragments.bed > ${name}_OCR_fragments.bed
  chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes`
  macs2 callpeak -f BEDPE -t ${name}_OCR_fragments.bed  --outdir ../3_peak -n ${name} -g ${chromsize} 2>&1 >>/dev/null | tee ../3_peak/${name}_MACS.out
  ShiftPairEnd.sh ${name}_OCR_fragments.bed ${fragment_length}
  n=`wc -l ${name}_OCR_fragments_shift.bed | cut -f 1 -d " "` && \
  c=`bc -l <<< "1000000 / $n"` && \
  genomeCoverageBed -bga -scale $c -i ${name}_OCR_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_OCR_fragments_shift.bdg && \
  bdg2bw.sh ${name}_OCR_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name} && \
  rm ${name}_OCR_fragments_shift.bdg ${name}_OCR_fragments_shift.bed
  cd ..
}

function SNP {
  # ----- bam split -----
  #cd 1_mapping
  #python3 ../${MY_PATH}/bam_split_snp.py ${SNP_info} ${name}.bam ${name}_C57BL.bam ${name}_DBA.bam > ${name}_split_snp.log
  #cd ..
  #
  #bamToBed -bedpe -i 1_mapping/${name1}.bam | awk '{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}' | sort -S 1% -k1,1 -k2,2n | uniq > 2_signal/${name1}_fragments.bed && \
  #awk '{print $3-$2}' 2_signal/${name1}_fragments.bed | sort -S 1% | uniq -c | sort -S 1% -k2,2g | awk 'BEGIN{print "fragment_length\tnumber"} {print $2"\t"$1}' > 2_signal/${name1}_fragments_length.txt &
  #bamToBed -bedpe -i 1_mapping/${name2}.bam | awk '{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}' | sort -S 1% -k1,1 -k2,2n | uniq > 2_signal/${name2}_fragments.bed && \
  #awk '{print $3-$2}' 2_signal/${name2}_fragments.bed | sort -S 1% | uniq -c | sort -S 1% -k2,2g | awk 'BEGIN{print "fragment_length\tnumber"} {print $2"\t"$1}' > 2_signal/${name2}_fragments_length.txt &
  #wait

  ## ----- peak_calling_pileup -----
  #chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
  #macs2 callpeak -f BEDPE -t 2_signal/${name1}_fragments.bed --outdir 2_peak -n ${name1} -g ${chromsize} -q 0.05 2>&1 >>/dev/null | tee 2_peak/MACS_{name1}.log && \
  #cd 2_signal && \
  fragment_length=`awk 'BEGIN{s=0;c=0} NR>1{s+=$1*$2;c+=$2} END{printf "%f", s/c}' ${name1}_fragments_length.txt` && \
  ShiftPairEnd.sh ${name1}_fragments.bed ${fragment_length} && \
  n=`wc -l ${name1}_fragments_shift.bed | cut -f 1 -d " "` && \
  c=`bc -l <<< "1000000 / $n"` && \
  genomeCoverageBed -bga -scale $c -i ${name1}_fragments_shift.bed -g ~/source/chrom_sizes/${genomeVersion}.chrom.sizes > ${name1}_fragments_shift.bdg && \
  bdg2bw.sh ${name1}_fragments_shift.bdg ~/source/chrom_sizes/${genomeVersion}_main.chrom.sizes ${name1} && \
  rm ${name1}_fragments_shift.bdg ${name1}_fragments_shift.bed && \
  cd .. &

  #chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
  #macs2 callpeak -f BEDPE -t 2_signal/${name2}_fragments.bed --outdir 2_peak -n ${name2} -g ${chromsize} -q 0.05 2>&1 >>/dev/null | tee 2_peak/MACS_${name2}.log && \
  #cd 2_signal && \
  fragment_length=`awk 'BEGIN{s=0;c=0} NR>1{s+=$1*$2;c+=$2} END{printf "%f", s/c}' ${name2}_fragments_length.txt` && \
  ShiftPairEnd.sh ${name2}_fragments.bed ${fragment_length} && \
  n=`wc -l ${name2}_fragments_shift.bed | cut -f 1 -d " "` && \
  c=`bc -l <<< "1000000 / $n"` && \
  genomeCoverageBed -bga -scale $c -i ${name2}_fragments_shift.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name2}_fragments_shift.bdg && \
  bdg2bw.sh ${name2}_fragments_shift.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name2} && \
  rm ${name2}_fragments_shift.bdg ${name2}_fragments_shift.bed && \
  cd .. &
  wait
  
  cd ..
}

# ----- running -----

mapping_filtering
short_fragments
#if [[ $# -eq 8 ]]; then
#  SNP
#fi
