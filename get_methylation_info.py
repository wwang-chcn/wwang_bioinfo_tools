#! /usr/bin/env python3

import os
import sys
from collections import defaultdict

from numpy import nanmean

self_name = sys.argv[0]
USAGE = f'{self_name} <name> <tssFile> <methylationFile> <span> <window> <step>'
if len(sys.argv) < 7:
    sys.stderr.write(f'No enough parameters!\n{USAGE}\n')
    sys.exit(1)

name            =     sys.argv[1]
tssFile         =     sys.argv[2]
methylationFile =     sys.argv[3]
span            = int(sys.argv[4])
window          = int(sys.argv[5])
step            = int(sys.argv[6])


def getMethylationInfo(tss, methylation, span, window, step):
    chrom, pos = tss
    coverage = []
    CpGnumber = []
    methylation_level = []
    unmethylatedCpG = []
    methylatedCpG = []
    mediumCpG = []
    for i in range(int(2*(span-window/2)/step+1)):
        local_methylation = []
        for j in range(window+1):
            if pos-span+i*step+j-1 in methylation[chrom]:
                local_methylation.append(methylation[chrom][pos-span+i*step+j-1])
        if len(local_methylation) == 0:
            coverage.append(float("nan"))
            CpGnumber.append(0)
            methylation_level.append(float("nan"))
            unmethylatedCpG.append(float("nan"))
            methylatedCpG.append(float("nan"))
            mediumCpG.append(float("nan"))
        else:
            coverage.append(nanmean([m[1] for m in local_methylation]))
            CpGnumber.append(len(local_methylation))
            methylation_level.append(1.0*sum([m[2] for m in local_methylation])/sum([m[1] for m in local_methylation]))
            unmethylatedCpG.append(sum(1 for m in local_methylation if m[0]<=0.2))
            methylatedCpG.append(sum(1 for m in local_methylation if m[0]>0.8))
            mediumCpG.append(sum(1 for m in local_methylation if m[0]<=0.8 and m[0]>0.2))

    return (methylation_level, unmethylatedCpG, coverage, CpGnumber, methylatedCpG, mediumCpG)



methylation = defaultdict(dict)
with open(methylationFile) as methylationFhd:
    for line in methylationFhd:
        if line[0] == "#" or not line:
            continue
        line = line.strip().split()
        methylation[line[0]][int(line[1])] = (float(line[3]), int(line[4]), int(line[5]))
    methylationFhd.close()


with open(                             tssFile ) as     tssFhd, \
     open(    name+"_methylationLevel.txt","w" ) as outputFhd1, \
     open(     name+"_unmethylatedCpG.txt","w" ) as outputFhd2, \
     open(            name+"_coverage.txt","w" ) as outputFhd3, \
     open(       name+"_methylatedCpG.txt","w" ) as outputFhd4, \
     open( name+"_middleMethylatedCpG.txt","w" ) as outputFhd5:
    for line in tssFhd:
        line = line.strip().split()
        (methylation_level, unmethylatedCpG, coverage, CpGnumber, methylatedCpG, mediumCpG) = getMethylationInfo((line[0], int(line[1])), methylation, span, window, step)
        outputFhd1.write("\t".join([str(m) for m in methylation_level])+"\n")
        outputFhd2.write("\t".join([str(m) for m in   unmethylatedCpG])+"\n")
        outputFhd3.write("\t".join([str(m) for m in          coverage])+"\n")
        outputFhd4.write("\t".join([str(m) for m in     methylatedCpG])+"\n")
        outputFhd5.write("\t".join([str(m) for m in         mediumCpG])+"\n")

    

