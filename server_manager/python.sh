#! /bin/bash

# Jan-22-2019

sudo apt update
sudo apt install python-setuptools python-pip python-numpy python-scipy python-pandas
sudo apt install default-libmysqlclient-dev # require by mysql-python
sudo pip install cutadapt # required by trim_galore
sudo pip install Cython # required by MACS2
sudo pip install pysam # required by HiC-Pro
sudo pip install hmmlearn # dependency: scikit_learn-0.20.2
sudo pip install beautifulsoup
sudo pip install bx-python # required by HiC-Pro
sudo pip install RSeQC
sudo pip install twobitreader
sudo pip install pyfasta
sudo pip install HTSeq
sudo pip install w3lib
sudo pip install scrapy
sudo pip install Biopython
sudo pip install jupyter
sudo pip install -U ggplot
sudo apt install python-virtualenv

wget https://github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz && tar zxvf MACS-1.4.2-1.tar.gz && cd MACS-1.4.2/ && sudo python setup.py install && cd ..
git clone https://github.com/taoliu/MACS.git && cd MACS/ && sudo python setup_w_cython.py install && cd ..
wget http://liulab.dfci.harvard.edu/CEAS/src/CEAS-Package-1.0.2.tar.gz && tar zxvf CEAS-Package-1.0.2.tar.gz && cd CEAS-Package-1.0.2/ && sudo python setup.py install && cd ..
wget http://cistrome.org/BETA/src/BETA_1.0.7.zip && unzip BETA_1.0.7.zip && cd CEAS-Package-1.0.2/ && sudo python setup.py install && cd ..
git clone https://github.com/hutuqiu/bseqc.git && cd bseqc && sudo python setup.py install && cd ..
hg clone https://bitbucket.org/cistrome/cistrome-applications-harvard && cd cistrome-applications-harvard/cistrome-extra-apps && sudo python setup.py install && cd ../mdseqpos && sudo cp lib/settings.py.example lib/settings.py && sudo python setup.py install && cd ../../
git clone https://github.com/ChengchenZhao/DrSeq2.git && cd DrSeq2 && sudo python setup.py install && cd ..
