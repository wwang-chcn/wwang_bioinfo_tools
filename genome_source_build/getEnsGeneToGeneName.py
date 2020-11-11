#! /usr/bin/env python

# Jun-28-2019
# Update: Nov-11-2020

import os
import sys
from collections import defaultdict

USAGE = f'{os.path.basename(sys.argv[0])} <genomeVersion>'

if len(sys.argv) < 2:
    sys.stdout.write('No enought arguments!\n')
    sys.stdout.write(USAGE+'\n')
    sys.exit(1)

genomeVersion = sys.argv[1]


with open(f'{genomeVersion}.ensemblToGeneName.txt') as fhd:
    ensToGeneSymbol = {line.strip().split()[0]: line.strip().split()[1] for line in fhd}

ensGeneToGeneSymbol = defaultdict(list)
with open(f'{genomeVersion}.ensGene.genePredExt') as fhd:
    for line in fhd:
        line = line.strip().split()
        if line[0] in ensToGeneSymbol:
            ensGeneToGeneSymbol[line[11]].append(ensToGeneSymbol[line[0]])
        else:
            ensGeneToGeneSymbol[line[11]].append('NA')

with open(f'{genomeVersion}.ensGeneToGeneSymbol.txt', 'w') as fhd:
    for gene in ensGeneToGeneSymbol:
        fhd.write(f'{gene}\t{",".join(set(ensGeneToGeneSymbol[gene]))}\n')