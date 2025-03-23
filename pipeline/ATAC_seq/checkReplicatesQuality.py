#!/usr/bin/env python3

# ------------------------------------
# Load Modules
# ------------------------------------
import os
import sys
import subprocess
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Helvetica'
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

# ------------------------------------
# Constant
# ------------------------------------
source_dir = os.path.expanduser('~/source/')


# ------------------------------------
# Sub Functions
# ------------------------------------
def sample_gen(running_file):
    """
    generate information about samples
    """
    with open(running_file) as fhd:
        for line in fhd:
            if not line.startswith('bash'):
                continue
            _, _, name, processer, genome_version, reads1_files, reads2_files, *_ = line.strip(
            ).split()
            yield name, genome_version, reads1_files, reads2_files


def gen_replicates_dict(sample_gen):
    """
    handle generated infomation about samples, output replicates (dict) and genome_version information
    """
    replicates = defaultdict(list)
    for name, genome_version, *_ in sample_gen:
        label, replicate = name.rsplit('_', 1)
        replicates[label].append(name)
    return replicates, genome_version


def check_target_regions_file(genome_version,
                              genome_bin_file=None,
                              tss_3000_file=None):
    """
    check whether target regions file exits, otherwise generate it
    """
    if genome_bin_file is None:
        genome_bin_file = f'{source_dir}/bySpecies/{genome_version}/{genome_version}.genome.5kb.bed'
    if tss_3000_file is None:
        tss_3000_file = f'{source_dir}/bySpecies/{genome_version}/{genome_version}.refGene.tss_3000.bed'

    # check target region file
    if not os.path.isfile(genome_bin_file):
        genome_size_file = f'{source_dir}/bySpecies/{genome_version}/{genome_version}_main.chrom.sizes'
        with open(genome_size_file) as input_fhd, \
             open(genome_bin_file, 'w') as output_fhd:
            index = 0
            for line in input_fhd:
                chrom, length = line.strip().split()
                for i in range(int(length) // 5000 + 1):
                    index += 1
                    output_fhd.write(
                        f'{chrom}\t{i*5000}\t{i*5000+5000}\tbin_{index}\n')
    if not os.path.isfile(tss_3000_file):
        def promoter_gen(fhd,promoter_range=3000):
            for index, line in enumerate(fhd):
                line = line.strip().split()
                name, chrom, strand, TSS, TTS = line[:5]
                if '_' in chrom or 'NR' in name:
                    continue
                TSS = int(TSS) if strand == '+' else int(TTS)
                new_line = [chrom, f'{max(TSS-promoter_range, 0):d}', f'{TSS+promoter_range:d}', f'{name}_{index+1}', '0', strand]
                yield '\t'.join(new_line) + '\n'
        refgene_file = f'{source_dir}/bySpecies/{genome_version}/{genome_version}.refGene.genePredExt'
        with open(refgene_file) as input_fhd, \
             open(tss_3000_file, 'w') as output_fhd:
            for line in promoter_gen(input_fhd):
                output_fhd.write(line)


def get_bigwig_mean(bigwig_file):
    """
    get the mean about a bigwig file
    """
    import subprocess
    fold = subprocess.check_output(f'bigWigInfo {bigwig_file} | grep mean',
                                   shell=True).decode()
    fold = float(fold.split()[1])
    return fold


def get_bigwig_files(samples):
    """
    get bigwig files for each sample
    """
    # under replicates_quality_check directory
    bigwig_files = []
    for sample in samples:
        bigwig_file = f'../2_signal/{sample}_uniq_SE_reads.bw'
        if os.path.isfile(bigwig_file):
            if not os.path.isfile(f'{sample}.bw'):
                os.symlink(src=bigwig_file, dst=f'{sample}.bw')
            bigwig_files.append(f'{sample}.bw')
    return bigwig_files


def avg_bed_signale_capture(name, target_regions_file, bigwig_files, labels, parallel_num=8):
    """
    perform bigWigAverageOverBed for given target_regions_file and bigwi_files
    """
    cmd = f'cut -f 1-4 {target_regions_file} > captures_regions.bed'
    subprocess.call(cmd, shell=True)
    cmd = ''
    for index, (label, bigwig_file) in enumerate(zip(labels, bigwig_files)):
        bw_scan_cmd = f'bigWigAverageOverBed {bigwig_file} captures_regions.bed {label}_signal.tsv &\n'
        cmd += bw_scan_cmd
        if (index + 1) % 8 == 0:
            cmd += 'wait\n'
    cmd += 'wait\n'
    with open('tmp.sh', 'w') as fhd:
        fhd.write(cmd)
    subprocess.call('bash tmp.sh'.split())


def process_signal(name, target_regions_file, bigwig_files, labels):
    """
    process the captured signal
    """
    capture_regions = pd.read_csv(target_regions_file, sep='\t', header=None)
    capture_signal = pd.DataFrame(index=capture_regions[3].to_list())
    for label, bigwig_file in zip(labels, bigwig_files):
        avg = get_bigwig_mean(bigwig_file)
        capture_signal[label] = pd.read_csv(
            f'{label}_signal.tsv',
            sep='\t',
            header=None,
            index_col=0,
            names=['size', 'covered', 'sum', 'mean0', 'mean'])['mean0'] / avg
    corr = capture_signal.corr()
    print(f'Correlation is:\n{corr}')
    corr_3 = {}
    mean_corr = lambda x: (x.values.sum() - 3) / 6
    for replicates_3 in combinations(corr.index, 3):
        sub_corr = corr.loc[replicates_3, replicates_3]
        corr_3[replicates_3] = mean_corr(sub_corr)
    corr_3 = sorted([[key, value] for key, value in corr_3.items()],
                    key=lambda x: x[1],
                    reverse=True)
    replicates, corr_value = corr_3[0]
    print(f'Replicates with highest correlation are: {replicates}')
    print(f'The mean correlation is: {corr_value}')


def main():
    os.makedirs('replicates_quality_check', exist_ok=True)
    replicates, genome_version = gen_replicates_dict(sample_gen('runned.sh'))
    # target regions files
    genome_bin_file = f'{source_dir}/bySpecies/{genome_version}/{genome_version}.genome.5kb.bed'
    tss_3000_file = f'{source_dir}/bySpecies/{genome_version}/{genome_version}.refGene.tss_3000.bed'
    check_target_regions_file(genome_version, genome_bin_file, tss_3000_file)
    # processing summary
    MY_PATH = os.path.dirname(__file__)
    cmd = f'python3 {MY_PATH}/ATAC_summary.py replicates_check'
    subprocess.call(cmd.split())
    os.rename(
        src='ATAC_summary_replicates_check_fragments_length.csv',
        dst=
        'replicates_quality_check/ATAC_summary_replicates_check_fragments_length.csv'
    )
    os.rename(
        src='ATAC_summary_replicates_check_mapping_info.csv',
        dst=
        'replicates_quality_check/ATAC_summary_replicates_check_mapping_info.csv'
    )
    #
    os.chdir('replicates_quality_check')
    for condition, samples in replicates.items():
        bigwig_files = get_bigwig_files(samples)
        avg_bed_signale_capture(condition, genome_bin_file, bigwig_files,
                                samples)
        process_signal(condition, genome_bin_file, bigwig_files, samples)
        avg_bed_signale_capture(condition, tss_3000_file, bigwig_files,
                                samples)
        process_signal(condition, tss_3000_file, bigwig_files, samples)
    os.chdir('../')


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write("User interrupts me! ;-) See you ^.^!")
        sys.exit(0)