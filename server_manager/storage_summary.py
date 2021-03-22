#! /usr/bin/env python3

import os
import pandas as pd
from functools import reduce
from pprint import pprint

storage = []
out = os.popen('du -h -d 1').read()

storage.append(pd.DataFrame({'user':[x[2:] for x in out.split()[1:-2:2]],'total':out.split()[0:-3:2]}))

file_regex_collection = {'sam': '.*\.sam', 'fastq': '.*\.f\(ast\)?q', 'wig': '.*\.wig', 'bdg': '.*\.b\(e\)?dg\(raph\)?', 'public_raw_data': '.*[SE]RR[0-9]+\(_[12]\)?\.f\(ast\)?q.gz'}
for file_regex_label in file_regex_collection:
    file_regex = file_regex_collection[file_regex_label]
    out  = os.popen(f'for i in *; do echo $i; find ./$i -type f -iregex "{file_regex}" -print0 | du -ch --files0-from=- | tail -n 1 ; done').read()
    storage.append(pd.DataFrame({'user':out.split('\n')[0:-1:2], file_regex_label: [x.split()[0] for x in out.split('\n')[1::2]]}))

storage_use = reduce(lambda left, right: pd.merge(left,right,left_on='user',right_on='user'),storage)
pprint(storage_use)

#find . -iregex ".*\.f\(ast\)?q" > fastq.files.txt
#find . -iregex ".*\.sam" > sam.files.txt
#find . -iregex ".*\.b\(e\)?dg\(raph\)?" > bdg.files.txt
#find . -iregex ".*[SE]RR[0-9]+\(_[12]\)?\.f\(ast\)?q.gz" > public_raw_data.files.txt