#! /bin/bash

# Jan-22-2019

# ----- server manage -----
sudo apt install traceroute
sudo apt install bc
sudo apt install mercurial # hg
sudo apt install unzip # require by others
sudo apt install autoconf # require by others
sudo apt install make # require by others
sudo apt install cmake # require by others
sudo apt install clang # require by others
sudo apt install gcc # require by others
sudo apt install g++ # require by others
sudo apt install libz-dev # require by bedtools
sudo apt install libbz2-dev # require by htslib
sudo apt install liblzma-dev # require by htslib
sudo apt install libncurses5-dev # require by samtools
sudo apt install libgsl-dev # require by bcftools, gfold
sudo apt install libboost-all-dev # require by tophat
sudo apt install libseqan2-dev # require by tophat
sudo apt install libxml2 libxml2-dev
sudo apt install texlive-latex-base # require by CAM
sudo apt install pigz # parallel gzip
sudo apt install exfat-fuse exfat-utils # mount exfat hardware

# ----- downloading -----

wget https://download.asperasoft.com/download/sw/connect/3.9.6/ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz && tar zxvf ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz && sudo bash ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh && chmod 755 ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh && sudo cp ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh /usr/local/bin

# ----- tools -----
# UCSC
mkdir UCSC_tools && cd UCSC_tools/ && rsync -aP rsync://hgdownload.cse.ucsc.edu/genome/admin/exe/linux.x86_64/ ./ && for i in *; do if [[ -x $i ]]; then sudo mv $i /usr/local/bin/ ; fi; done && sudo mv blat/blat blat/gfClient blat/gfServer /usr/local/bin/ && cd ..
# htslib, samtools, bcftools
git clone https://github.com/samtools/htslib.git && cd htslib/ && autoheader && autoconf && ./configure && make && sudo make install && cd ..
git clone https://github.com/samtools/samtools.git && cd samtools/ && autoheader && autoconf -Wno-syntax && ./configure && make && sudo make install && cd ..
git clone https://github.com/samtools/bcftools.git && cd bcftools/ && autoheader && autoconf -Wno-syntax && ./configure && make && sudo make install && cd ..
# bedtools require python, htslib
git clone https://github.com/arq5x/bedtools2.git && cd bedtools2/ && make && sudo make install && cd ..

# ncbi software
mkdir ncbi && cd ncbi/
git clone https://github.com/ncbi/ncbi-vdb.git && cd ncbi-vdb/ && ./configure && make && sudo make install && cd .. # require ngs-sdk, require by ngs-python & need re-login
git clone https://github.com/ncbi/ngs.git && cd ngs/ && ./configure && make -C ngs-sdk && make -C ngs-java && make -C ngs-python && sudo make -C ngs-sdk install && sudo make -C ngs-java install && sudo make -C ngs-python install && cd .. # require java, python, python-setuptools; required by hisat2,
git clone https://github.com/ncbi/sra-tools.git && cd sra-tools && ./configure && make && sudo make install && sudo mv /usr/local/ncbi/sra-tools/bin/* /usr/local/bin/ && cd ..
git clone https://github.com/ncbi/ngs-tools.git && cd ngs-tools/ && ./configure && 
cd ..



# ----- jhu softwares -----
# bowtie
git clone https://github.com/BenLangmead/bowtie && cd bowtie/ && make NO_TBB=1 && sudo make install && cd ..
# bowtie2 
git clone https://github.com/BenLangmead/bowtie2.git && cd bowtie2/ && make NO_TBB=1 && sudo make install && cd ..
# tophat
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && tar zxvf tophat-2.1.1.Linux_x86_64.tar.gz && cd tophat-2.1.1.Linux_x86_64/ && sudo cp bam2fastx bam_merge bed_to_juncs contig_to_chr_coords fix_map_ordering gtf_juncs gtf_to_fasta juncs_db long_spanning_reads map2gtf prep_reads sam_juncs samtools_0.1.18 segment_juncs sra_to_solid tophat tophat-fusion-post tophat2 tophat_reports /usr/local/bin/ && cd ..
# cufflinks
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && tar zxvf cufflinks-2.2.1.Linux_x86_64.tar.gz && cd cufflinks-2.2.1.Linux_x86_64 && sudo cp cuffcompare  cuffdiff  cufflinks  cuffmerge  cuffnorm  cuffquant  gffread  gtf_to_sam /usr/local/bin/ && cd ..
# hisat2
wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-source.zip && unzip hisat2-2.1.0-source.zip && cd hisat2-2.1.0/ && make USE_SRA=1 NCBI_NGS_DIR=/usr/local/ngs/ngs-sdk NCBI_VDB_DIR=/usr/local/ncbi/ncbi-vdb && sudo cp hisat2 hisat2-align-s hisat2-align-l hisat2-build hisat2-build-s hisat2-build-l hisat2-inspect hisat2-inspect-s hisat2-inspect-l /usr/local/bin/ && cd ..
# stringtie
git clone https://github.com/gpertea/stringtie && cd stringtie/ && make release && sudo cp prepDE.py stringtie /usr/local/bin/ && cd ..

# ----- other mapping softwares -----
# bwa
git clone https://github.com/lh3/bwa && cd bwa && make && sudo cp bwa /usr/local/bin/ && cd ..
# STAR
git clone https://github.com/alexdobin/STAR.git && cd STAR/source && make STAR && cd ../../ && sudo cp STAR/bin/Linux_x86_64/* /usr/local/bin/

# ----- DNA methylation softwares -----
# bismark
git clone https://github.com/FelixKrueger/Bismark.git && cd Bismark/ && sudo cp NOMe_filtering bam2nuc bismark bismark2bedGraph bismark2report bismark2summary bismark_genome_preparation bismark_methylation_extractor copy_bismark_files_for_release.pl coverage2cytosine deduplicate_bismark filter_non_conversion /usr/local/bin/ && cd ..
# moabs
git clone https://github.com/sunnyisgalaxy/moabs.git && cd moabs/bin/ && sudo cp bsmap mcall mcomp moabs numCI preprocess_novoalign.sh redepth.pl sam2bam.sh /usr/local/bin/ && cd ../..

# ----- meme suite -----
# ghost
wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs927/ghostscript-9.27-linux-x86_64.tgz && tar zxvf ghostscript-9.27-linux-x86_64.tgz && sudo cp ghostscript-9.27-linux-x86_64/gs-927-linux-x86_64 /usr/local/bin/gs
wget http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz && tar zxvf meme-5.0.5.tar.gz && cd meme-5.0.5/ && ./configure --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt && make && make test && sudo make install && cd .. && sudo cp /home/compbio/bin/ame /home/compbio/bin/centrimo /home/compbio/bin/dreme /home/compbio/bin/dreme-py3 /home/compbio/bin/dust /home/compbio/bin/fimo /home/compbio/bin/glam2 /home/compbio/bin/glam2scan /home/compbio/bin/gomo /home/compbio/bin/mast /home/compbio/bin/mcast /home/compbio/bin/meme /home/compbio/bin/meme-chip /home/compbio/bin/momo /home/compbio/bin/purge /home/compbio/bin/spamo /home/compbio/bin/tomtom ./scripts/ama-qvalues ./scripts/beeml2meme ./scripts/chen2meme ./scripts/dreme_xml_to_html ./scripts/dreme_xml_to_txt ./scripts/elm2meme ./scripts/fasta-center ./scripts/fasta-dinucleotide-shuffle ./scripts/fasta-dinucleotide-shuffle-py3 ./scripts/fasta-dinucleotide-shuffle-py3.py ./scripts/fasta-dinucleotide-shuffle.py ./scripts/fasta-fetch ./scripts/fasta-grep ./scripts/fasta-hamming-enrich ./scripts/fasta-hamming-enrich-py3 ./scripts/fasta-hamming-enrich-py3.py ./scripts/fasta-hamming-enrich.py ./scripts/fasta-make-index ./scripts/fasta-most ./scripts/fasta-re-match ./scripts/fasta-subsample ./scripts/fasta-unique-names ./scripts/glam2html ./scripts/glam2psfm ./scripts/glam2scan2html ./scripts/hart2meme-bkg ./scripts/hartemink2psp ./scripts/iupac2meme ./scripts/jaspar2meme ./scripts/mast_xml_to_html ./scripts/mast_xml_to_txt ./scripts/matrix2meme ./scripts/meme-chip_html_to_tsv ./scripts/meme-rename ./scripts/meme_xml_to_html ./scripts/nmica2meme ./scripts/plot-usage ./scripts/plotgen ./scripts/priority2meme ./scripts/psp-gen ./scripts/rna2meme ./scripts/rsat-retrieve-seq ./scripts/rsat-supported-organisms ./scripts/scpd2meme ./scripts/sd ./scripts/sites2meme ./scripts/taipale2meme ./scripts/tamo2meme ./scripts/tomtom_xml_to_html ./scripts/transfac2meme ./scripts/uniprobe2meme ./scripts/usage-reports /usr/local/bin/

# ----- gfold -----
sudo hg clone https://bitbucket.org/feeldead/gfold && cd gfold && sudo make && sudo cp gfold /usr/local/bin/ && cd ..

# ----- HiC-pro -----
# require R, bx-python, pysam, ggplot2
git clone https://github.com/nservant/HiC-Pro.git && cd HiC-Pro/
# edit config-install.txt install position
# edit ./scripts/install/install_dependencies.sh samtools version check
# samver=`samtools 2>&1 | grep Version | cut -d" " -f2 | cut -d"-" -f1`
make configure && sudo make install && cd ..

# ----- Docker -----
sudo apt install docker.io
sudo systemctl start docker
sudo systemctl enable docker

# ----- Virtual display software ------
sudo apt install xvfb
