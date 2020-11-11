#! /usr/bin/env python

import os
import sys

USAGE = f'{os.path.basename(sys.argv[0])} <genomeVersion>'

if len(sys.argv) < 2:
    sys.stdout.write('No enought arguments!\n')
    sys.stdout.write(USAGE+'\n')
    sys.exit(1)

genomeVersion = sys.argv[1]

with open(f'{genomeVersion}.fa') as inputFhd, \
     open(f'{genomeVersion}_main.fa', 'w') as outputFhd:
    flag = False
    for line in inputFhd:
        if line[0] == '>':
            if '_' in line:
                flag = False
            else:
                flag = True
                outputFhd.write(line)
        else:
            if flag:
                outputFhd.write(line)