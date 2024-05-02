#!/usr/bin/env python3

# ------------------------------
# Import Modules
# ------------------------------
import os
import subprocess
import sys
from optparse import OptionParser

import numpy as np
from scipy.ndimage import gaussian_filter


# ------------------------------
# Sub Funcitons
# ------------------------------
def prepare_optparser():
    '''\
    Prepare optparser object. New options will be added in thisfunction first.
    '''
    script_name = os.path.basename(sys.argv[0])
    usage = f'USAGE: {script_name} <-n name> <-p peak_files+> <-b bigwig_files+>'
    description = 'merge_peak -- merge peaks from replicates'

    # option processor
    optparser = OptionParser(version=f'{script_name} 0.1', description=description, usage=usage, add_help_option=False)

    # basic setting
    optparser.add_option('-h', '--help', action='help', help='Show this help message and exit.')
    optparser.add_option('-n', '--name',dest='name', type='string',\
                         help='Name for this sample.')
    optparser.add_option('-p', '--peak',dest='peak', type='string', action='append',\
                         help='Peaks files for each replicates.')
    optparser.add_option('-q', '--q-cutoff',dest='q_cutoff', type='float',\
                         help='Qvalue cutoff.',default='6')
    optparser.add_option('-f', '--f-cutoff',dest='f_cutoff', type='float',\
                         help='Fold cutoff.',default='6')
    optparser.add_option('-b', '--bigwig',dest='bigwig', type='string', action='append',\
                         help='Bigwig files for each replicates. Must be same number as peak files.')
    optparser.add_option('-g', '--genomesize',dest='gs', type='string',\
                         help='Genomesize file.')

    return optparser


def opt_validate(optparser):
    '''Validate options from a OptParser object.
    Ret: Validated options object.
    '''
    (options, args) = optparser.parse_args()

    # input name must be given
    if not options.name:
        sys.stdout.write('Input name be given!\n')
        optparser.print_help()
        sys.exit(1)

    # input bigwig and bed files must be given
    if not (len(options.peak) > 1 and len(options.bigwig) > 1):
        sys.stdout.write('At least 2 replicate peaks and bigWig files must be given!\n')
        optparser.print_help()
        sys.exit(1)
    if not len(options.peak) == len(options.bigwig):
        sys.stdout.write('The number of peaks files and bigWig files must be same!\n')
        optparser.print_help()
        sys.exit(1)

    # check each input file
    for i in range(len(options.peak)):
        if not os.path.isfile(options.peak[i]):
            sys.stdout.write(f'Can not find the peaks file: {options.peak[i]}.\n')
            optparser.print_help()
            sys.exit(1)

    for i in range(len(options.bigwig)):
        if not os.path.isfile(options.bigwig[i]):
            sys.stdout.write(f'Can not find the bigwig file: {options.bigwig[i]}.\n')
            optparser.print_help()
            sys.exit(1)

    return options


def peak_intersect(options):
    cmd = f'cat {" ".join(options.peak)} | awk "{{if(\$7>{options.f_cutoff} && \$9>{options.q_cutoff}) print}}" | cut -f 1-3 | grep -v chrM | sort -k1,1 -k2,2n | mergeBed -i - > temp_{options.name}_merged_peaks.bed'
    subprocess.call(cmd, shell=True)
    cmd = f'cat temp_{options.name}_merged_peaks.bed | intersectBed -u -a - -b {" | intersectBed -u -a - -b ".join(options.peak)} > temp_{options.name}_intersect_peaks.bed'
    subprocess.call(cmd, shell=True)
    os.remove(f'temp_{options.name}_merged_peaks.bed')


def merge_bw(options):
    cmd = f'bigWigMerge -threshold=-1 {" ".join(options.bigwig)} temp1_{options.name}_merged.bdg'
    subprocess.call(cmd.split())
    cmd = f'sort -k1,1 -k2,2n temp1_{options.name}_merged.bdg > temp2_{options.name}_merged.bdg'
    subprocess.call(cmd, shell=True)
    cmd = f'bedGraphToBigWig temp2_{options.name}_merged.bdg {options.gs} {options.name}.bw'
    subprocess.call(cmd.split())
    os.remove(f'temp1_{options.name}_merged.bdg')
    os.remove(f'temp2_{options.name}_merged.bdg')


def gen_single_peak_signal(peak_file, signal_file):
    with open(peak_file) as fhd:
        for line in fhd:
            chrom, start, end, *_ = line.strip().split('\t')
            start, end = int(start), int(end)
            cmd = f'bigWigSummary {signal_file} {chrom} {start} {end} {end-start}'
            output = [float(x) for x in subprocess.check_output(cmd.split()).decode().strip().split('\t')]
            yield chrom, start, end, gaussian_filter(output, sigma=(end - start) * 0.05, mode='constant')


def find_summits(options):
    with open(f'{options.name}_peaks.bed', 'w') as fhd1, \
         open(f'{options.name}_summits.bed', 'w') as fhd2:
        for index, (chrom, start, end, signal) in enumerate(
                gen_single_peak_signal(f'temp_{options.name}_intersect_peaks.bed', f'{options.name}.bw')):
            fhd1.write(f'{chrom}\t{start}\t{end}\t{options.name}_peak_{index+1}\t{np.mean(signal):.5f}\t.\n')
            fhd2.write(
                f'{chrom}\t{start+np.argmax(signal)}\t{start+np.argmax(signal)+1}\t{options.name}_peak_{index+1}\t{np.mean(signal):.5f}\t.\n'
            )
    os.remove(f'temp_{options.name}_intersect_peaks.bed')
    os.remove(f'{options.name}.bw')


# ------------------------------
# Main Funcitons
# ------------------------------
def main():

    # read the options and validate them
    options = opt_validate(prepare_optparser())

    peak_intersect(options)
    merge_bw(options)
    find_summits(options)


# ------------------------------
# Program Running
# ------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you!\n')
        sys.exit(0)
