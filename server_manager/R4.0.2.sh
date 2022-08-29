#! /bin/bash

# Aug-4-2020

# required for compile
sudo apt install tcl8.6-dev tk8.6-dev # tcl/tk
sudo apt install texinfo texlive texlive-fonts-extra

wget https://cran.r-project.org/src/base/R-4/R-4.0.2.tar.gz
tar zxvf R-4.0.2.tar.gz
./configure --prefix=/usr/lib/R/4.0.2 --enable-R-shlib --enable-memory-profiling
make
sudo make install

cd /usr/bin/
sudo ln -s /usr/lib/R/4.0.2/bin/R R4.0.2 && sudo ln -s /usr/lib/R/4.0.2/bin/Rscript Rscript4.0.2


# check png, pdf et al.
#capabilities()
