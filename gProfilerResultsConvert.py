#! /usr/bin/env python

import os, sys
from itertools import dropwhile
import csv

USAGE = '{} <outputFile> <gProfilerResults> <idToGeneSymbolFile>'.format(os.path.basename(sys.argv[0]))

if len(sys.argv) < 4:
    sys.stdout.write('No enought parameters, exit!\n')
    sys.stdout.write(USAGE+'\n')
    sys.exit(1)

outputFile = sys.argv[1]
gProfilerResultsFile = sys.argv[2]
idToGeneSymbolFile = sys.argv[3]

with open(idToGeneSymbolFile) as fhd:
    idToGeneSymbol = {line.strip().split()[0]: line.strip().split()[1] for line in fhd}

def idConvert(rows):
    for row in rows:
        row[-1] = ','.join(idToGeneSymbol[geneID] if geneID in idToGeneSymbol else 'NA' for geneID in row[-1].split(','))
        term_size, query_size, intersection_size, effective_domain_size = int(row[5]), int(row[6]), int(row[7]), int(row[8])
        row.insert(5,(1.0*intersection_size/query_size)/(1.0*term_size/effective_domain_size))
        yield row

with open(gProfilerResultsFile) as inputFhd, \
     open(outputFile, 'w') as outputFhd:
    input_csv = csv.reader(dropwhile(lambda x: x.startswith('#'),inputFhd))
    output_csv = csv.writer(outputFhd)
    header = next(input_csv)
    header.insert(5,'Enrichment')
    output_csv.writerow(header)
    output_csv.writerows(idConvert(input_csv))

