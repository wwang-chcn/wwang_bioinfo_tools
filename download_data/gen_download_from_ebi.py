#! /usr/bin/env python3

import os
import sys
import csv

USAGE = '{} <proj_file> <output_file>'.format(os.path.basename(sys.argv[0]))

proj_file = sys.argv[1]
output_file = sys.argv[2]


def get_fastq_file_path(row):
    fastq_file_paths, fastq_files = [], []
    if 'fastq_aspera' in row:
        for download_link in row.get('fastq_aspera').split(';'):
            fastq_file_paths.append(download_link.split(':')[1])
            fastq_files.append(download_link.rsplit('/',1)[1])
    elif 'fastq_ftp' in row:
        for download_link in row.get('fastq_ftp').split(';'):
            fastq_file_paths.append(download_link.split('/')[1])
            fastq_files.append(download_link.rsplit('/',1)[1])
    else:
        return None, None
    
    return fastq_file_paths, fastq_files


def gen_download_script(rows):
    for row in rows:
        fastq_file_paths, fastq_files = get_fastq_file_path(row)
        if not fastq_file_paths:
            line_text = r'\t'.join(row.values())
            sys.stderr.write(f'Unexpected line: "{line_text}", skip.\n')
            continue
        for fastq_file_path, fastq_file in zip(fastq_file_paths, fastq_files):
            yield f'if [ ! -f {fastq_file} ]; then ~/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:{fastq_file_path} .; fi'
            yield f'if [ ! -f {fastq_file} ]; then wget ftp://ftp.sra.ebi.ac.uk/{fastq_file_path}; fi'


with open(proj_file) as input_fhd, \
     open(output_file, 'w') as output_fhd:
    intput_csv = csv.DictReader(input_fhd, delimiter='\t')
    output_fhd.write('\n'.join(gen_download_script(intput_csv))+'\n')
