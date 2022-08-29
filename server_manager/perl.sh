#! /bin/bash

# Jan-23-2019

sudo cpan Module::Build

curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.0.tar.gz -o trim_galore.tar.gz && tar xvzf trim_galore.tar.gz && sudo ln -s /home/compbio/ServerPackages/TrimGalore-0.6.0/trim_galore /usr/local/bin/
mkdir homer && cd homer && wget http://homer.ucsd.edu/homer/configureHomer.pl && perl configureHomer.pl -install homer && sudo ln -s /home/compbio/ServerPackages/homer/bin/* /usr/local/bin/ && cd ../

