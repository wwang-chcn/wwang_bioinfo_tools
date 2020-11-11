#! /usr/bin/env python

# Jun-28-2019

import os, sys
from collections import defaultdict

USAGE = '{} <genomeVersion>'.format(os.path.basename(sys.argv[0]))

if len(sys.argv) < 2:
    sys.stdout.write('No enought arguments!\n')
    sys.stdout.write(USAGE+'\n')
    sys.exit(1)

genomeVersion = sys.argv[1]


with open('{}.ensemblToGeneName.txt'.format(genomeVersion)) as fhd:
    ensToGeneSymbol = {line.strip().split()[0]: line.strip().split()[1] for line in fhd}

ensGeneToGeneSymbol = defaultdict(list)
with open('{}.ensGene.genePredExt'.format(genomeVersion)) as fhd:
    for line in fhd:
        line = line.strip().split()
        if line[0] in ensToGeneSymbol:
            ensGeneToGeneSymbol[line[11]].append(ensToGeneSymbol[line[0]])
        else:
            ensGeneToGeneSymbol[line[11]].append('NA')

with open('{}.ensGeneToGeneSymbol.txt'.format(genomeVersion), 'w') as fhd:
    for gene in ensGeneToGeneSymbol:
        fhd.write('{}\t{}\n'.format(gene,','.join(set(ensGeneToGeneSymbol[gene]))))