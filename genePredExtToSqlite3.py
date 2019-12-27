#! /usr/bin/env python

import os, sys
import sqlite3

USAGE = '{} <genePredExt> <outputDatabase> [wigFile]'.format(os.path.basename(sys.argv[0]))

if len(sys.argv) < 3:
    sys.stdout.write('No enought arguments!\n')
    sys.stdout.write(USAGE+'\n')
    sys.exit(1)

genePredExtFile = sys.argv[1]
sqlite3File = sys.argv[2]

conn = sqlite3.connect(sqlite3File)
c = conn.cursor()
c.execute("CREATE TABLE GeneTable (chrom,name,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,name2)")

with open(genePredExtFile) as fhd:
    for line in fhd:
        line = line.strip().split()
        exonstart = line[8].split(',')
        exonend = line[9].split(',')
        newline = [line[1], line[0], line[2], int(line[3]), int(line[4]), int(line[5]), int(line[6]), int(line[7]), ', '.join(exonstart).strip(), ', '.join(exonend).strip(), line[11]]
        c.execute("insert into GeneTable (chrom,name,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts, exonEnds, name2) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (newline))

conn.commit()
conn.close()

if len(sys.argv) < 4:
    sys.stdout.write('WIG file not given, bgCorrected database will not generate.\n')
    sys.exit(0)

import subprocess
p = subprocess.Popen('build_genomeBG -g {} -w {} -o {}.bgCorrected.sq3'.format(sqlite3File,sys.argv[3],sqlite3File.rsplit('.',1)[0]), shell=True)
sts = os.waitpid(p.pid, 0)


