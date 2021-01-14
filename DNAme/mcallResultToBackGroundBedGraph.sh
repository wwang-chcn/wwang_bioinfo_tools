#! /bin/bash

function print_help {
    echo "USAGE: $0 <mcall.out> [threshold=0]"
}

if [[ -e $1 ]]
then
    if [[ -n $2 ]] # has threshold value
    then
        if [[ $2 =~ ^[1-9][0-9]?$ ]] # input is a positive interge(right input)
        then
            threshold=$2
        else
            print_help; exit 2
        fi
    else
        threshold=1
    fi
else
    print_help; exit 1
fi

awk -v threshold=$threshold '$1 !~/#/{if($5>=threshold && NF==15) print $1"\t"$2"\t"$3"\t1"}' $1 | sort -k1,1 -k2,2n