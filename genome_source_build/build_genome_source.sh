#! /bin/bash

# Update: Nov-11-2020

USAGE="<genomeVersion>"

MY_PATH="`dirname \"$0\"`"
genomeVersion=${1}

rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/bigZips/${genomeVersion}.2bit .
twoBitToFa ${genomeVersion}.2bit ${genomeVersion}.fa

python ${MY_PATH}/getMainFa.py ${genomeVersion}

# mapping index
bowtie2-build --threads 16 ${genomeVersion}.fa ${genomeVersion}
hisat2-build -p 16 ${genomeVersion}.fa ${genomeVersion}
bwa index ${genomeVersion}.fa

bowtie2-build --threads 16 ${genomeVersion}_main.fa ${genomeVersion}_main
hisat2-build -p 16 ${genomeVersion}_main.fa ${genomeVersion}_main
bwa index ${genomeVersion}_main.fa

mkdir -p bismark_index/{raw,main}
cd bismark_index/raw && ln -s ../../${genomeVersion}.fa .
cd ../main && ln -s ../../${genomeVersion}_main.fa .
cd ../..
bismark_genome_preparation --parallel 8 --single_fasta ./bismark_index/raw
bismark_genome_preparation --parallel 8 --single_fasta ./bismark_index/main


rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/bigZips/${genomeVersion}.chrom.sizes .

# annotation
# refGene
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/refGene.txt.gz ${genomeVersion}.refGene.txt.gz
gunzip ${genomeVersion}.refGene.txt.gz
cat ${genomeVersion}.refGene.txt | cut -f 2-11 | genePredToGtf -utr file stdin ${genomeVersion}.refGene.gtf
python ${MY_PATH}/correctGtfGeneID.py ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.gtf
cut -f 2-11 ${genomeVersion}.refGene.txt > ${genomeVersion}.refGene.genePred
cut -f 2-16 ${genomeVersion}.refGene.txt > ${genomeVersion}.refGene.genePredExt
genePredToBed ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.bed


# ensGene
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/ensGene.txt.gz ${genomeVersion}.ensGene.txt.gz
if [[ -e ${genomeVersion}.ensGene.txt.gz ]]; then
    gunzip ${genomeVersion}.ensGene.txt.gz
    cat ${genomeVersion}.ensGene.txt | cut -f 2-11 | genePredToGtf -utr file stdin ${genomeVersion}.ensGene.gtf
    python ${MY_PATH}/correctGtfGeneID.py ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.gtf
    cut -f 2-11 ${genomeVersion}.ensGene.txt > ${genomeVersion}.ensGene.genePred
    cut -f 2-16 ${genomeVersion}.ensGene.txt > ${genomeVersion}.ensGene.genePredExt
    genePredToBed ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.bed
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/ensemblToGeneName.txt.gz ${genomeVersion}.ensemblToGeneName.txt.gz
    gunzip ${genomeVersion}.ensemblToGeneName.txt.gz
    python ${MY_PATH}/getEnsGeneToGeneName.py ${genomeVersion}
fi

python ${MY_PATH}/genePredExtToSqlite3.py ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.sq3
if [[ -e ${genomeVersion}.ensGene.txt.gz ]]; then
    python ${MY_PATH}/genePredExtToSqlite3.py ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.sq3
fi

