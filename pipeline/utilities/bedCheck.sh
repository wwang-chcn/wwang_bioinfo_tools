#!/usr/bin/env bash

bedFile=${1}
bbFile=${bedFile::(${#bedFile}-2)}b

#check if bedFile exist
if [[ ! -e ${bedFile} && -e ${bbFile} ]]; then
    bigBedToBed ${bbFile} ${bedFile}