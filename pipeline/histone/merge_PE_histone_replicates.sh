#! /bin/bash

# May-28-2021

function print_help {
    echo "USAGE: $0 <name> <genomeVersion> <replicates,+> <narrow/broad> [downSampling]"
    echo "name: the sample name"
    echo "genomeVersion: the genome version processed for these samples"
    echo "replicates: comma separated list of replicates"
    echo "narrow/broad: peak type"
    echo "downSampling: target sampling fragments for each sample. optional"
}

function bedToBigWig {
    n=`wc -l ${name}_fragments.bed | cut -f 1 -d " "`
    c=`bc -l <<< "1000000 / $n"`
    genomeCoverageBed -bga -scale $c -i ${name}_fragments.bed -g ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes > ${name}_fragments.bdg
    ${MY_PATH}/../utilities/bdg2bw.sh ${name}_fragments.bdg ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${name}
    rm ${name}_fragments.bdg
}

function compress_bed {
    bedFile=${1}
    genomeVersion=${2}
    
    if [[ ! -e ${bedFile::(${#bedFile}-2)}b ]]; then
        col=`head -1 ${bedFile} | awk '{print NF}'`
        plus=`bc <<< "$col -3"`
        intersectBed -a ${bedFile} -b <(awk '{print $1"\t0\t"$2}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes) -wa -f 1.00 | sort -k1,1 -k2,2n > ${bedFile}.tmp
        bedToBigBed -type=bed3+${plus} ${bedFile}.tmp ~/source/bySpecies/${genomeVersion}/${genomeVersion}.chrom.sizes ${bedFile::(${#bedFile}-2)}b
        rm ${bedFile} ${bedFile}.tmp
    fi
}

function peak_calling {
    pc_name=${1}
    genomeVersion=${2}
    peakType=${3}
    
    cd 5_merged_sample
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes`
    if [[ ! -e ${pc_name}_fragments.bed ]]; then
        bigBedToBed ${pc_name}_fragments.bb ${pc_name}_fragments.bed
    fi
    chromsize=`awk 'BEGIN{s=0} {s+=$2} END{print s}' ~/source/bySpecies/${genomeVersion}/${genomeVersion}_main.chrom.sizes` && \
    if [[ ${peakType} == "narrow" ]]; then
        macs3 callpeak -f BEDPE -t ${pc_name}_fragments.bed -n ${pc_name} -g ${chromsize} --nomodel --shift 37 --extsize 73 --keep-dup all 2>&1 >>/dev/null | tee MACS_${pc_name}.out
    elif [[ ${peakType} == "broad" ]]; then
        macs3 callpeak -f BEDPE -t ${pc_name}_fragments.bed -n ${pc_name} -g ${chromsize} --nomodel --shift 37 --extsize 73 --broad --keep-dup all 2>&1 >>/dev/null | tee MACS_${pc_name}.out
    else
        echo "peak type should be narrow or broad, exit!"
        exit 1
    fi
    cd ..
}

function clearning_up {
    for i in ${fragments_array[@]}; do
        if [[ -e 2_signal/${i}_fragments.bed ]]; then
            rm 2_signal/${i}_fragments.bed
        fi
    done
    if [[ -e 5_merged_sample/${name}_fragments.bed ]]; then
        rm 5_merged_sample/${name}_fragments.bed
    fi
}

if [[ $# -lt 3 ]]; then
    echo "No enought parameters!"
    print_help
    exit 1
fi

MY_PATH="`readlink -f $(dirname \"$0\")`"

mkdir -p 5_merged_sample/

name=${1}
genomeVersion=${2}

# get fragment files
IFS=',' read -r -a reads_file_array <<< ${3}
for i in ${reads_file_array[@]}; do
    if [[ ! -e 2_signal/${i}_fragments.bed ]]; then
        if [[ ! -e 2_signal/${i}_fragments.bb ]]; then
            echo "fragments file for sample ${i} do not exist, exit!"
            exit 1
        fi
        bigBedToBed 2_signal/${i}_fragments.bb 2_signal/${i}_fragments.bed
    fi
    fragments_array+=("2_signal/${i}_fragments.bed")
done
peakType=${4}


# merge fragment files
cat ${fragments_array[@]} | sort -k1,1 -k2,2n > 5_merged_sample/${name}_fragments.bed.tmp
cd 5_merged_sample/
if [[ $# -eq 5 ]]; then
    mv ${name}_fragments.bed.tmp ${name}_fragments.bed
else
    n=`wc -l ${name}_fragments.bed.tmp | cut -f 1 -d " "`
    ratio=`bc -l <<< "${5} / $n"`
    awk -v ratio=${ratio} 'BEGIN{srand();}{if(rand() <= ratio) print $0}' ${name}_fragments.bed.tmp > ${name}_fragments.bed
    rm ${name}_fragments.bed.tmp
fi
bedToBigWig
cd ..

peak_calling ${name} ${genomeVersion} ${peakType}
compress_bed 5_merged_sample/${name}_fragments.bed ${genomeVersion}

clearning_up
