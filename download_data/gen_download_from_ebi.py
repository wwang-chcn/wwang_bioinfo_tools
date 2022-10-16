#! /usr/bin/env python3

import os
import sys
import csv

USAGE = '{} <proj_file> <output_file>'.format(os.path.basename(sys.argv[0]))

proj_file = sys.argv[1]
output_file = sys.argv[2]

def gen_download_script(rows):
    for row in rows:
        for fastq_file in row.get('fastq_ftp', None).split(';'):
            fastq_file_path = fastq_file.split('/',1)[1]
            fastq_file = fastq_file.rsplit('/',1)[1]
            yield f'if [ ! -f {fastq_file} ]; then ~/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:{fastq_file_path} .; fi'
            yield f'if [ ! -f {fastq_file} ]; then wget ftp://ftp.sra.ebi.ac.uk/{fastq_file_path}; fi'


with open(proj_file) as input_fhd, \
     open(output_file, 'w') as output_fhd:
    intput_csv = csv.DictReader(input_fhd, delimiter='\t')
    output_fhd.write('\n'.join(gen_download_script(intput_csv))+'\n')
