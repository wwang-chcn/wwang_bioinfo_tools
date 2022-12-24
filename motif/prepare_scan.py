#! /usr/bin/env python

import os
import sys

self_name = os.path.basename(sys.argv[0])
USAGE = f'{self_name} <motif_database> <genome>'

if len(sys.argv) < 3:
    sys.stdout.write(f'No enough parameters!\n{USAGE}\n')

motif_database, genome = sys.argv[1], sys.argv[2]

CMD = ""
output = os.popen(f'egrep -i ^motif {motif_database}').read()
n = 0
for line in output.split("\n"):
    if not line:
        continue
    n += 1
    line = line.split()
    motifID, motifName = line[1], line[2].replace("(", "_").replace(")", "")
    CMD += """fimo --oc raw_results/{0}_{1} --max-stored-scores 1000000 --motif {1} {2} {3} && awk 'NR>1 && /^[^#]/ && NF>9{{print $3"\\t"$4-1"\\t"$5"\\t"$10"\\t"$8"\\t"$6}}' raw_results/{0}_{1}/fimo.tsv | sort -k1,1 -k2,2n > results/{0}_{1}.bed &\n""".format(
        motifName, motifID, motif_database, genome)
    if n % 20 == 0:
        CMD += "wait\n"

with open('motifScan.sh', 'w') as fhd:
        fhd.write(f'''#! /bin/bash

{CMD}
wait
''')


CMD1 = ""
n = 0
for line in output.split("\n"):
    if not line:
        continue
    n += 1
    line = line.split()
    motifID, motifName = line[1], line[2].replace("(", "_").replace(")", "")
    CMD1 += """awk 'NR>1 && /^[^#]/ && NF>9{{print $3"\\t"$4-1"\\t"$5"\\t"$10"\\t"$8"\\t"$6}}' raw_results/{0}_{1}/fimo.tsv | sort -k1,1 -k2,2n > results/{0}_{1}.bed &\n""".format(
        motifName, motifID, motif_database, genome)
    if n % 20 == 0:
        CMD1 += "wait\n"
with open('motifScan_temp.sh', 'w') as fhd:
        fhd.write(f'''#! /bin/bash

{CMD1}
wait
''')
