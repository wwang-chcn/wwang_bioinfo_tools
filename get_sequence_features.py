#! /usr/bin/env python3

# Jun-6-2021

import os
import sys
import csv
from optparse import OptionParser
import numpy as np

# ------------------------------------
# Sub Functions
# ------------------------------------
def prepare_optparser():
    
    '''Prepare optparser object. New options will be added in this
    function first.
    '''
    
    program_name = os.path.basename(sys.argv[0])
    usage = 'usage: %prog <-g genome_file> <-b bed_file> <-o output_name> <-w window>'
    description = 'Get sequence feature for given regions.'
    
    # option processor
    optparser = OptionParser(version='%prog 0.1', description=description, usage=usage, add_help_option=False)
    optparser.add_option('-h','--help', action='help', help='Show this help message and exit.')
    optparser.add_option('-g','--genome', dest='genome', type='string',\
                         help='Required. Genome sequence file. Can be 2bit or fasta format, auto detect by suffix.')
    optparser.add_option('-b','--bed',dest='bed',type='string',\
                         help='Required. Bed file for regions to be capture.')
    optparser.add_option('-o','--output',dest='output',type='string',\
                         help='Required. Output file basename.')
    optparser.add_option('-w','--window',dest='window',type='int',\
                         help='Required. Sliding window size.')
    return optparser


def opt_validate(optparser):
    
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    
    (options, args) = optparser.parse_args()
    
    # input genome sequence file must be given
    if not options.genome:
        sys.stdout.write('Error: Input genome sequence file must be given!\nExit\n')
        optparser.print_help()
        sys.exit(1)
        if not os.path.isfile(os.path.expanduser(options.genome)):
            sys.stdout.write(f'Error: Input genome sequence file: {options.genome} can not be found!\nExit\n')
            optparser.print_help()
            sys.exit(1)
        suffix = options.genome.split('.')[-1]
        if not suffix in ['2bit','fa','fasta']:
            sys.stdout.write('Error: Input genome sequence file must be 2bit or fasta format!\nExit\n')
            optparser.print_help()
            sys.exit(1)
        if suffix == '2bit':
            sys.stdout.write('Info: Input genome sequence file in 2bit format.\n')
        else:
            sys.stdout.write('Info: Input genome sequence file in fasta format.\n')
    
    # input bed file must be given
    if not options.bed:
        sys.stdout.write('Error: Input bed file must be given!\nExit\n')
        optparser.print_help()
        sys.exit(1)
        if not os.path.isfile(os.path.expanduser(options.bed)):
            sys.stdout.write(f'Error: Input bed file: {options.bed} can not be found!\nExit\n')
            optparser.print_help()
            sys.exit(1)
    
    # output name must be given
    if not options.output:
        sys.stdout.write('Error: Output name must be given!\nExit\n')
        optparser.print_help()
        sys.exit(1)
    
    # window and step must be set
    if not options.window:
        sys.stdout.write(f'Error: window is not set.\n')
        optparser.print_help()
        sys.exit(1)
    
    return options


def cal_sequence_feature(sequence):
    '''
    Input: upper cased sequence
    Output: GC content, CpG density, CpG ratio
    '''
    n = len(sequence)
    if n == 0:
        return np.nan, np.nan, np.nan
    G = sequence.count('G')
    C = sequence.count('C')
    CpG = sequence.count('CG')
    GC, CpG_density, CpG_ratio = 1.0 * (G + C) / n, 1.0 * CpG / n, 1.0 * CpG * n / C / G if C and G else np.nan
    return GC, CpG_density, CpG_ratio


# ------------------------------------
# Main function
# ------------------------------------
def main():
    
    # read the options
    options = opt_validate(prepare_optparser())
    
    # load sequence file
    suffix = options.genome.split('.')[-1]
    if suffix == '2bit':
        from twobitreader import TwoBitFile
        genome = TwoBitFile(os.path.expanduser(options.genome))
    else:
        from pyfasta import Fasta
        genome = Fasta(os.path.expanduser(options.genome))

    # calculate matrix size
    n = 1
    with open(os.path.expanduser(options.bed)) as input_fhd:
        line = next(input_fhd).strip().split()
        start, end = int(line[1]), int(line[2])
        for line in input_fhd:
            n += 1
    points = (end - start) // options.window
    GC_matrix = np.zeros(shape=(n,points))
    CpG_density_matrix = np.zeros(shape=(n,points))
    CpG_ratio_matrix = np.zeros(shape=(n,points))
    
    # get sequence features matrix
    i = 0
    with open(os.path.expanduser(options.bed)) as input_fhd:
        for line in input_fhd:
            if not line:
                continue
            line = line.strip().split()
            chrom, start, end, name = line[:4]
            start, end = int(start), int(end)
            start = max(start, 0)
            sequence = genome[chrom][start:end].upper()
            for j in range((end-start)//options.window):
                GC, CpG_density, CpG_ratio = cal_sequence_feature(sequence[j*options.window:j*options.window+options.window])
                GC_matrix[i,j] = GC
                CpG_density_matrix[i,j] = CpG_density
                CpG_ratio_matrix[i,j] = CpG_ratio
            i += 1

    # output
    np.savetxt(f'{os.path.expanduser(options.output)}_GC_content.csv', GC_matrix, delimiter=',', fmt='%.3f')
    np.savetxt(f'{os.path.expanduser(options.output)}_CpG_density.csv', CpG_density_matrix, delimiter=',', fmt='%.3f')
    np.savetxt(f'{os.path.expanduser(options.output)}_CpG_ratio.csv', CpG_ratio_matrix, delimiter=',', fmt='%.3f')

# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you ^.^!\n')
        sys.exit(0)

