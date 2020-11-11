#! /usr/bin/env python

# Jun-28-2019
# Update: Nov-11-2020
# More pythonic

import os
import sys
from collections import defaultdict

USAGE = f'{os.path.basename(sys.argv[1])} <genePredExt> <gtf>'

if len(sys.argv) < 3:
    sys.stdout.write('No enought arguments!\n')
    sys.stdout.write(USAGE+'\n')
    sys.exit(1)

genePredExtFile = sys.argv[1]
gtfFile = sys.argv[2]

with open(genePredExtFile) as fhd:
    gene = {line.strip().split()[0]: line.strip().split()[11] for line in fhd}

with open(gtfFile) as inputFhd, \
     open('gtf.tmp', 'w') as outputFhd:
    for line in inputFhd:
        line = line.strip().split('\t')
        attribute = line[-1].split()
        try:
            attribute[1] = '"'+gene[attribute[1][1:-2]]+'";'
        except KeyError:
        	transcript = attribute[1][1:-2].rsplit('_',1)
        	attribute[1] = '"'+gene[transcript[0]]+'_'+transcript[1]+'";'
        outputFhd.write('{}\t{}\n'.format('\t'.join(line[:-1]),' '.join(attribute)))

os.rename('gtf.tmp',gtfFile)