#! /usr/bin/env python3

# Jan-14-2020

import os
import sys
from optparse import OptionParser
import numpy as np
from collections import defaultdict


# ------------------------------------
# Sub Functions
# ------------------------------------
def prepare_optparser():
    
    '''Prepare optparser object. New options will be added in this
    function first.
    '''
    
    program_name = os.path.basename(sys.argv[0])
    usage = 'usage: %prog <-m meth_file> <-b bed_file> [-o output_name] [-c coverage]'
    description = 'get_methylation_level.'
    
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

    # output name
    if not options.output:
        meth_name = options.meth.split('.')[0]
        options.output = [f'{bed_file[:-3]}_{meth_name}_methylation_level.txt' for bed_file in options.bed]
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
            methylation[line[0]][int(line[1])] = (float(line[3]), int(line[4]), int(line[5]))
    return methylation


def get_region_meth(region, methylation, coverage):
    
    (chrom, start, end) = region
    m, total, n = 0, 0, 0
    for pos in range(start-1,end):
        if pos in methylation[chrom]:
            n += 1
            m += methylation[chrom][pos][2]
            total += methylation[chrom][pos][1]
    return np.nan if n == 0 or total / n < coverage else m / total


def get_file_meth(bed_file, output_file, methylation, coverage):
    basename = bed_file.rsplit('.',1)[0]
    with open(bed_file) as intput_fhd, \
         open(output_file, 'w') as output_fhd:
        line_n = 0
        for line in intput_fhd:
            if not line:
                continue
            line_n += 1
            line = line.strip().split()
            chrom, start, end = line[:3]
            start, end = int(start), int(end)
            name = line[3] if len(line) > 3 else f'R{line_n:d}'
            value = get_region_meth((chrom, start, end), methylation, coverage)
            output_fhd.write(f'{chrom}\t{start:d}\t{end:d}\t{name}\t{value:.3f}\t.\n')
 

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
        get_file_meth(bed_file, output_file, methylation, options.coverage)


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you ^.^!\n')
        sys.exit(0)
