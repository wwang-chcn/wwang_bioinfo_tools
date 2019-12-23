#! /usr/bin/env python3

import os
import sys
import pandas as pd

USAGE = '{} <proj_file> <output_file>'.format(os.path.basename(sys.argv[0]))

proj_file = sys.argv[1]
output_file = sys.argv[2]

def gen_download_script(rows):
    for row in rows:
        for fastq_file in row.fastq_ftp.split(';'):
            fastq_file_path = fastq_file.split('/',1)[1]
            yield f'~/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:{fastq_file_path} .'


proj = pd.read_csv(proj_file, sep='\t')

with open(output_file, 'w') as fhd:
    fhd.write('\n'.join(gen_download_script(proj.itertuples())))
