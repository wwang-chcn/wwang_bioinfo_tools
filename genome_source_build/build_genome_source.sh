#! /bin/bash

# Update: Nov-11-2020

USAGE="bash $0 <genomeVersion>"

function print_help {
    echo "USAGE: $USAGE"
}

if [[ $# -lt 1 ]]; then
    echo No enought parameters!
    print_help
    exit 1
fi

MY_PATH="`readlink -f $(dirname \"$0\")`"
genomeVersion=${1}

if [[ ! -e ${genomeVersion}.2bit ]]; then
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/bigZips/${genomeVersion}.2bit .
fi
twoBitToFa ${genomeVersion}.2bit ${genomeVersion}.fa

python ${MY_PATH}/getMainFa.py ${genomeVersion}

# mapping index
bowtie2-build --threads 16 ${genomeVersion}.fa ${genomeVersion}
hisat2-build -p 16 ${genomeVersion}.fa ${genomeVersion}
bwa index ${genomeVersion}.fa

bowtie2-build --threads 16 ${genomeVersion}_main.fa ${genomeVersion}_main
hisat2-build -p 16 ${genomeVersion}_main.fa ${genomeVersion}_main
bwa index ${genomeVersion}_main.fa

twoBitInfo ${genomeVersion}.2bit ${genomeVersion}.chrom.sizes
grep -v "_" ${genomeVersion}.chrom.sizes > ${genomeVersion}_main.chrom.sizes


# annotation
# refGene
if [[ ! -e ${genomeVersion}.refGene.txt.gz ]]; then
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/refGene.txt.gz ${genomeVersion}.refGene.txt.gz
fi
gunzip ${genomeVersion}.refGene.txt.gz
cut -f 2-11 ${genomeVersion}.refGene.txt > ${genomeVersion}.refGene.genePred
cut -f 2-16 ${genomeVersion}.refGene.txt > ${genomeVersion}.refGene.genePredExt
genePredToGtf -utr file ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.gtf
genePredToBed ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.bed


# ensGene
if [[ ! -e ${genomeVersion}.ensGene.txt.gz ]]; then
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/ensGene.txt.gz ${genomeVersion}.ensGene.txt.gz
fi
if [[ -s ${genomeVersion}.ensGene.txt.gz ]]; then
    gunzip ${genomeVersion}.ensGene.txt.gz
    cut -f 2-11 ${genomeVersion}.ensGene.txt > ${genomeVersion}.ensGene.genePred
    cut -f 2-16 ${genomeVersion}.ensGene.txt > ${genomeVersion}.ensGene.genePredExt
    genePredToGtf -utr file ${genomeVersion}.ensGene.genePred ${genomeVersion}.ensGene.gtf
    python ${MY_PATH}/correctGtfGeneID.py ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.gtf
    genePredToBed ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.bed
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/ensemblToGeneName.txt.gz ${genomeVersion}.ensemblToGeneName.txt.gz
    gunzip ${genomeVersion}.ensemblToGeneName.txt.gz
    python ${MY_PATH}/getEnsGeneToGeneName.py ${genomeVersion}
fi

python ${MY_PATH}/genePredExtToSqlite3.py ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.sq3
if [[ -e ${genomeVersion}.ensGene.txt.gz ]]; then
    python ${MY_PATH}/genePredExtToSqlite3.py ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.sq3
fi

