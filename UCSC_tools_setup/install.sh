#! /bin/bash

function print_help {
    echo "$0 <dst_dir>"
    echo "Install the exectueable files into dst_dir."
}


dst_dir=$1
for file in * blat/*; do if [[ -x $file && ! -d $file ]]; then echo $file; fi; done > ucsc_tools_list.txt
for file in `cat ucsc_tools_list.txt`; do cp $file ${dst_dir}/; done