#! /usr/bin/env python3

import os
import subprocess
from random import choice
from collections import defaultdict

# load blast_db
blast_db_dir = '~/source/blast_db/'
for file_name in os.listdir(os.path.expanduser(blast_db_dir)):
    src = os.path.expanduser(blast_db_dir + file_name)
    dst = file_name
    if not os.path.isfile(dst):
        os.symlink(src, dst)

samples = [
    'Oblongsmart2-_BKDL202560375-1a-N701-N501_1.fq.gz'
]
select_reads_number = 10000
for index, sample in enumerate(samples):
    cmd = f'''zcat {sample} | head -{select_reads_number*4} | awk 'NR % 4 == 1{{print ">"substr($$0,2,length($$0))}} NR % 4 == 2{{print $$0}}' > test.fa'''
    subprocess.call(cmd,shell=True)
    cmd = f"blastn -db nt -query test.fa -out test_output{index}.txt -outfmt '6 qseqid evalue bitscore staxid ssciname scomname'"
    subprocess.call(cmd,shell=True)

def alignment_process(lines):
    max_score = max(float(line[2]) for line in lines)
    organisms = [line[-1] for line in lines if float(line[2]) == max_score]
    return choice(organisms)

# plotting
with sns.axes_style('whitegrid'), sns.plotting_context('paper'):
    for index, sample in enumerate(samples):
        with open(f'test_output{index}.txt') as fhd:
            spicies = defaultdict(int)
            line = next(fhd)
            line = line.strip().split('\t')
            lines, name = [line], line[0]
            for line in fhd:
                line = line.strip().split('\t')
                if line[0] != name:
                    spicies[alignment_process(lines)] += 1
                    lines, name = [line], line[0]
                else:
                    lines.append(line)
            spicies[alignment_process(lines)] += 1
        organisms = [x[0] for x in  sorted(spicies.items(), key=lambda x: x[1], reverse=True)[:10]]
        percentages = [x[1]/10000*100 for x in sorted(spicies.items(), key=lambda x: x[1], reverse=True)[:10]]
        spicies = {'organism': organisms, 'pecentage': percentages}
        fig, ax = plt.subplots()
        sns.barplot(x='organism',y='pecentage',data=spicies,ax=ax)
        title = sample.split('_')[1]
        ax.set_xticklabels(spicies['organism'],rotation='vertical')
        with open(f'../1_mapping/Mapping_oblong_20200302_rep{index+1}.log') as fhd:
            for line in fhd:
                if 'Overall alignment rate' in line:
                    ar = line.strip().split()[-1]
        ax.set_title(f'{title} (alignment rate: {ar})')
        ax.set_ylabel('Percentage')
        fig.tight_layout()
        fig.savefig(f'smartseq_sample{index+1}_species.pdf')
        plt.close()