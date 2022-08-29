#! /bin/bash

# Jan-23-2019

# ----- R basic -----
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
#sudo echo deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/ >> /etc/apt/sources.list
sudo apt update
sudo apt upgrade
sudo apt install r-base
sudo apt install r-base-dev
sudo apt install littler python-rpy2 jags
sudo apt-get install libgdal-dev libproj-dev # for r-java

# RStudio
sudo apt install gdebi-core
wget https://download2.rstudio.org/server/trusty/amd64/rstudio-server-1.2.1335-amd64.deb
sudo gdebi rstudio-server-1.2.1335-amd64.deb

# require for curl
sudo apt install libcurl4-openssl-dev
# require for ssl
sudo apt install libssl-dev
# require for units
sudo apt install libudunits2-dev
# require for XML, xml2, litter
sudo apt install libxml2-dev
# require for export
sudo apt install libfontconfig1-dev libcairo2-dev libglu1-mesa-dev
