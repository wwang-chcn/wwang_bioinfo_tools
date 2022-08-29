#! /bin/bash

# Jan-23-2019

sudo apt update
# java
sudo apt install default-jre
# javac
sudo apt install default-jdk
# fastqc required by trim_galore
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && unzip fastqc_v0.11.8.zip && chmod 755 FastQC/fastqc && sudo ln -s /home/compbio/ServerPackages/FastQC/fastqc /usr/local/bin/fastqc # softlink required
# picard
wget https://github.com/broadinstitute/picard/releases/download/2.18.24/picard.jar
# chromHMM
wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip && unzip ChromHMM.zip
