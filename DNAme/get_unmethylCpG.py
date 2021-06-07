#! /usr/bin/env python3

# Jun-6-2021

import os
import sys
import csv
from optparse import OptionParser
import numpy as np
from collections import defaultdict
from contextlib import contextmanager


# ------------------------------------
# Sub Functions
# ------------------------------------
def prepare_optparser():
    
    '''Prepare optparser object. New options will be added in this
    function first.
    '''
    
    program_name = os.path.basename(sys.argv[0])
    usage = 'usage: %prog <-m meth_file> <-b bed_file> [-o output_name] [-c coverage] [-w window] [-s step]'
    description = 'Get unmethylCpG number.'
    
    # option processor
    optparser = OptionParser(version='%prog 0.1', description=description, usage=usage, add_help_option=False)
    optparser.add_option('-h','--help', action='help', help='Show this help message and exit.')
    optparser.add_option('-m','--meth', dest='meth', type='string',\
                         help='Methylation info file.')
    optparser.add_option('-b','--bed',dest='bed',type='string', action='append',\
                         help='Bed file(s) for regions to be capture.')
    optparser.add_option('-o','--output',dest='output',type='string', action='append',\
                         help='Names for output file(s).')
    optparser.add_option('-c', '--coverage',dest='coverage',type='int',\
                         help='Coverage threshold for region to be consider.\nDefault is 5.', default=5)
    optparser.add_option('-w','--window',dest='window',type='int',\
                         help='Sliding window size. Use the whole region if not set.')
    optparser.add_option('-s','--step',dest='step',type='int',\
                         help='Step for sliding window. Required when -w is set. Ignored if -w not set.')
    return optparser


def opt_validate(optparser):
    
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    
    (options, args) = optparser.parse_args()
    
    # input methylation info file must be given
    if not options.meth:
        sys.stdout.write('Input methylation info file must be given!\n')
        optparser.print_help()
        sys.exit(1)
    if not os.path.isfile(options.meth):
        sys.stdout.write(f'Input methylation info file: {options.meth} can not be found!\n')
        optparser.print_help()
        sys.exit(1)
    
    # input bed file(s) must be given
    if not options.bed:
        sys.stdout.write('Input bed file(s) must be given!\n')
        optparser.print_help()
        sys.exit(1)
    for bed_file in options.bed:
        if not os.path.isfile(bed_file):
            sys.stdout.write(f'Input bed file: {bed_file} can not be found!\n')
            optparser.print_help()
            sys.exit(1)
    
    # step should be set when window is set
    if options.window and not options.step:
        sys.stdout.write(f'Warning: -w window is set but -s step not set, ignored.\n')
    if options.step and not options.window:
        sys.stdout.write(f'Warning: -s step is set but -w window not set, ignored.\n')
    
    # output name
    if not options.output:
        meth_name = options.meth.split('.')[0]
        options.output = [f'{bed_file[:-3]}_{meth_name}_unmethylCpG.csv' for bed_file in options.bed]
    elif len(options.output) != len(options.bed):
        sys.stdout.write('The number of input bed files and number of output files are not equal!\n')
        optparser.print_help()
        sys.exit(1)
    
    return options


def load_meth(meth_file):
    methylation = defaultdict(dict)
    # {chrom: pos: {(ratio,totalC,methC)}}
    with open(meth_file) as fhd:
        for line in fhd:
            if line[0] == '#' or not line:
                continue
            line = line.strip().split()
            if len(line) < 15: # handle extra line in the bottom of XX.bam.G.bed
                continue
            methylation[line[0]][int(line[1])] = (float(line[3]), int(line[4]), int(line[5]))
    return methylation


def _get_region_unmethylCpG(region, methylation, coverage):
    (chrom, start, end) = region
    um = 0
    for pos in range(start-1,end):
        if pos in methylation[chrom]:
            if methylation[chrom][pos][1] >=5 and methylation[chrom][pos][2] <= 0.2:
                um += 1
    return um


def get_region_unmethylCpG(region, methylation, coverage, options):
    '''
    '''
    (chrom, start, end) = region
    if options.window and options.step:
        region_um = []
        for i in range((end - start - options.window) // options.step):
            _start = start + i * options.step
            _end = _start + options.window
            region_um.append(_get_region_unmethylCpG((chrom, _start, _end), methylation, coverage))
        return np.nanmax(region_um) if region_um else np.nan
    else:
        return _get_region_unmethylCpG(region, methylation, coverage)


@contextmanager
def smart_write_open(file_name):
    fhd = gzip.open(file_name, 'wt') if file_name.endswith('.gz') else open(file_name, 'w')
    yield fhd
    fhd.close()


def get_file_meth(bed_file, output_file, methylation, coverage, options):
    basename = bed_file.rsplit('.',1)[0]
    with open(bed_file) as intput_fhd, \
         smart_write_open(output_file) as output_fhd:
        output_csv = csv.writer(output_fhd)
        output_csv.writerow(['chrom','start','end','name','value','strand'])
        line_n = 0
        for line in intput_fhd:
            if not line:
                continue
            line_n += 1
            line = line.strip().split()
            chrom, start, end = line[:3]
            start, end = int(start), int(end)
            name = line[3] if len(line) > 3 else f'R{line_n:d}'
            strand = line[5] if len(line) > 5 else '.'
            value = get_region_unmethylCpG((chrom, start, end), methylation, coverage, options)
            output_line = [chrom, start, end, name, f'{value:.3f}', strand]
            output_csv.writerow(output_line)
 

# ------------------------------------
# Main function
# ------------------------------------
def main():

    # read the options
    options = opt_validate(prepare_optparser())

    # load methylation information
    methylation = load_meth(options.meth)

    # get methylation information for each file
    for bed_file, output_file in zip(options.bed, options.output):
        get_file_meth(bed_file, output_file, methylation, options.coverage, options)


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you ^.^!\n')
        sys.exit(0)
