#!/usr/bin/env python3

# ------------------------------
# Import Modules
# ------------------------------
import os
import sys
from optparse import OptionParser
import subprocess
from multiprocessing import cpu_count, Pool
import time
from time import localtime, strftime
import numpy as np
import pandas as pd
import gzip
from shutil import rmtree


# ------------------------------
# Sub Funcitons
# ------------------------------
def prepare_optparser():
    '''\
    Prepare optparser object. New options will be added in thisfunction first.
    '''
    script_name = os.path.basename(sys.argv[0])
    usage = f'USAGE: {script_name} <-b bedFile> <-w bigWigFiles+> [-n name] [-p process] [-s datapoints] [-m model]'
    description = 'bigWigSignalCapture -- Capture signal form bigWigFiles.'
    
    # option processor
    optparser = OptionParser(version=f'{script_name} 0.1',
                             description=description,
                             usage=usage,
                             add_help_option=False)
    
    # basic setting
    optparser.add_option('-h',
                         '--help',
                         action='help',
                         help='Show this help message and exit.')
    optparser.add_option('-n', '--name',dest='name', type='string',\
                         help='Name of this run. If not given, the time will be used.')
    optparser.add_option('-p', '--process',dest='process', type='int',\
                         help='Number of subprocess. If not given, number of computer processor will be used.',default=cpu_count())
    optparser.add_option('-s', '--datapoints',dest='datapoints', type='int',\
                         help='The capture region will broken into dataPoints equal parts. Default is 1.', default=1)
    optparser.add_option('-b', '--bedFile',dest='bedFile', type='string',\
                         help='BedFile for capture.')
    optparser.add_option('-w', '--bigWig',dest='bigWigFiles', type='string', action='append',\
                         help='bigWigFiles for capture.')
    optparser.add_option('-m', '--model',dest='model', type='string',\
                         help='Capture model: security or speed. Default is Security',default='security')
    optparser.add_option('--strand',dest='strand', action='store_true',\
                         help='If Strand specific')
    
    return optparser


def opt_validate(optparser):
    '''Validate options from a OptParser object.
    Ret: Validated options object.
    '''
    (options, args) = optparser.parse_args()
    
    # input bigwig and bed files must be given
    if not (options.bedFile and options.bigWigFiles):
        sys.stdout.write('Input bigwig and bed files must be given!\n')
        optparser.print_help()
        sys.exit(1)
    
    # check each input file
    for i in range(len(options.bigWigFiles)):
        if not os.path.isfile(options.bigWigFiles[i]):
            sys.stdout.write(
                f'Can not find the bigwig file: {options.bigWigFiles[i]}.\n')
            optparser.print_help()
            sys.exit(1)
        options.bigWigFiles[i] = os.path.expanduser(options.bigWigFiles[i])
    
    if not os.path.isfile(options.bedFile):
        sys.stdout.write(f'Can not find the bigwig file: {options.bedFile}.\n')
        optparser.print_help()
        sys.exit(1)
    
    # check capture model
    if not (options.model == 'security' or options.model == 'speed'):
        sys.stdout.write('Invalid capture model! Must be security or speed.\n')
        optparser.print_help()
        sys.exit(1)
    
    # get name
    if not options.name:
        options.name = f'bigWigSignalCapture_{time.strftime("%Y.%b.%d.%H-%M-%S", time.localtime())}'
    
    return options


def prepare_capture(options):
    """load bed file and bigwig files into tmp dir"""
    # load bed
    os.makedirs(options.name)
    tmp_fhd = []
    for i in range(options.datapoints):
        tmp_fhd.append(open(f'{options.name}/tmp_{i}.bed', 'w'))
    with open(options.bedFile) as fhd:
        for line in fhd:
            chrom, start, end, name, value, strand, *_ = line.strip().split()
            start, end = int(start), int(end)
            boundaries = np.linspace(start, end, options.datapoints + 1).astype('int')
            for i in range(options.datapoints):
                tmp_fhd[i].write(f'{chrom}\t{boundaries[i]}\t{boundaries[i+1]}\t{name}\t{value}\t{strand}\n')
    for i in range(options.datapoints):
        tmp_fhd[i].close()
    # creat temp dirctory
    CMD = f'cd {options.name}\n'
    # load bigwig files
    for i in range(len(options.bigWigFiles)):
        if options.bigWigFiles[i].startswith('/'):
            CMD += f'ln -s {options.bigWigFiles[i]} .\n'
        else:
            CMD += f'ln -s ../{options.bigWigFiles[i]} .\n'
    os.system(CMD)


def avg_signal_over_bed(bed_file, bigwig_file, output_file):
    """bigWigAverageOverBed.
Parameters:
    bed_file: bed file name.
    bigwig_file: bigwig file name.
    output_file: output file name.
Return
    None."""
    import subprocess
    from distutils.spawn import find_executable
    bigWigAverageOverBed_path = f'{find_executable("bigWigAverageOverBed")}' if find_executable("bigWigAverageOverBed") else '../bigWigAverageOverBed'
    CMD = f'{bigWigAverageOverBed_path} {bigwig_file} {bed_file} {output_file}'
    # cmd_output = subprocess.check_output(CMD.split()).decode()
    sp = subprocess.Popen(CMD.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
    out, err = sp.communicate()


def capture_unit(file_names):
    """wrapper for parallel running"""
    avg_signal_over_bed(*file_names)


def reverse_by_flag(row, flag):
    return row.to_list()[::-1] if flag else row.to_list()


def capture_signal(options):
    print(strftime("Start capture: %a, %d %b %Y %H:%M:%S", localtime()))
    # get input & output for each basic capture unit
    file_names = []
    for bigwig_file_index, bigwig_file in enumerate(options.bigWigFiles):
        for datapoint_index in range(options.datapoints):
            bed_file = f'tmp_{datapoint_index}.bed'
            output_file = f'tmp_{bigwig_file_index}_{datapoint_index}.tsv'
            file_names.append([bed_file, bigwig_file, output_file])
    # capture using multiple processer
    os.chdir(options.name)
    with Pool(processes=options.process) as pool:
        pool.map(capture_unit, file_names)
    os.chdir('..')
    # load signal
    print(strftime("Start post-processing: %a, %d %b %Y %H:%M:%S", localtime()))
    with open(options.bedFile) as fhd:
        names = [line.strip().split()[3] for line in fhd]
        if options.strand:
            reversed_flag_dict = {'-': True, '+': False}
            reversed_flags = [reversed_flag_dict[line.strip().split()[5]] for line in fhd]
    for bigwig_file_index, bigwig_file in enumerate(options.bigWigFiles):
        signal_matrix = pd.DataFrame(index=names)
        for datapoint_index in range(options.datapoints):
            output_file = f'tmp_{bigwig_file_index}_{datapoint_index}.tsv'
            signal = pd.read_csv(f'{options.name}/{output_file}', header=None, sep='\t', names=['name','size','covered','sum','mean0','mean'])
            signal_matrix[datapoint_index] = signal['mean']
        if options.strand:
            arr = []
            for reversed_flag, (_, row) in zip(reversed_flags, signal_matrix.iterrows()):
                arr.append(reverse_by_flag(row, reversed_flag))
            signal_matrix = pd.DataFrame(data=arr, index=names, columns=np.arange(options.datapoints))
        signal_matrix.to_csv(f'signal_{options.name}_siteprof{bigwig_file_index+1}.csv.gz')
    rmtree(options.name)
    print(strftime("End post-processing: %a, %d %b %Y %H:%M:%S", localtime()))


# ------------------------------
# Main Funcitons
# ------------------------------
def main():
    
    # read the options and validate them
    options = opt_validate(prepare_optparser())
    
    # prepare capture
    prepare_capture(options)
    
    # capture signal
    capture_signal(options)


# ------------------------------
# Program Running
# ------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you!\n')
        sys.exit(0)
