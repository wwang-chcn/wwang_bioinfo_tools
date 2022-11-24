#! /usr/bin/env python3

# Jun-6-2021

import os
import sys
import csv
from optparse import OptionParser

# ------------------------------------
# Sub Functions
# ------------------------------------
def prepare_optparser():
    
    '''Prepare optparser object. New options will be added in this
    function first.
    '''
    
    program_name = os.path.basename(sys.argv[0])
    usage = 'usage: %prog <-g genome_file> <-b bed_file> <-o output_name> [-w window] [-s step]'
    description = 'Get sequence feature for given regions.'
    
    # option processor
    optparser = OptionParser(version='%prog 0.1', description=description, usage=usage, add_help_option=False)
    optparser.add_option('-h','--help', action='help', help='Show this help message and exit.')
    optparser.add_option('-g','--genome', dest='genome', type='string',\
                         help='Genome sequence file. Can be 2bit or fasta format, auto detect by suffix.')
    optparser.add_option('-b','--bed',dest='bed',type='string',\
                         help='Bed file for regions to be capture.')
    optparser.add_option('-o','--output',dest='output',type='string',\
                         help='Output file.')
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
    
    # output file must be given
    if not options.output:
        sys.stdout.write('Error: Output file must be given!\nExit\n')
        optparser.print_help()
        sys.exit(1)
    if os.path.isfile(os.path.expanduser(options.output)):
        sys.stdout.write(f'Warning: Output file: {options.output} alread exist!\n')
        optparser.print_help()
        sys.exit(1)
    
    # step should be set when window is set
    if options.window and not options.step:
        sys.stdout.write(f'Warning: -w window is set but -s step not set, ignored.\n')
    if options.step and not options.window:
        sys.stdout.write(f'Warning: -s step is set but -w window not set, ignored.\n')

    return options


def cal_sequence_feature(sequence, options):
    '''
    Input: upper cased sequence
    Output: GC content, CpG density, CpG ratio
    '''
    n = len(sequence)
    if options.window and options.step:
        GC, CpG_density, CpG_ratio = [], [], []
        for i in range((n - options.window) // options.step):
            sequence_ = sequence[options.step*i:options.step*i+options.window]
            G = sequence.count('G')
            C = sequence.count('C')
            CpG = sequence.count('CG')
            GC.append(1.0 * (G + C) / n)
            CpG_density.append(1.0 * CpG / n)
            CpG_ratio.append(1.0 * CpG * n / C / G)
        GC, CpG_density, CpG_ratio = max(GC), max(CpG_density), max(CpG_ratio)
    else:
        G = sequence.count('G')
        C = sequence.count('C')
        CpG = sequence.count('CG')
        GC, CpG_density, CpG_ratio = 1.0 * (G + C) / n, 1.0 * CpG / n, 1.0 * CpG * n / C / G
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
    with open(os.path.expanduser(options.bed)) as input_fhd, \
         open(os.path.expanduser(options.output), 'w') as output_fhd:
        output_csv = csv.writer(output_fhd)
        output_csv.writerow(['name','GC content','CpG density', 'CpG ratio'])
        for line in input_fhd:
            if not line:
                continue
            line = line.strip().split()
            chrom, start, end, name = line[:4]
            start, end = int(start), int(end)
            start = max(start, 0)
            sequence = genome[chrom][start:end].upper()
            GC, CpG_density, CpG_ratio = cal_sequence_feature(sequence, options)
            output_csv.writerow([name, f'{GC:.3f}', f'{CpG_density:.3f}', f'{CpG_ratio:.3f}'])


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write('User interrupts me! ;-) See you ^.^!\n')
        sys.exit(0)

