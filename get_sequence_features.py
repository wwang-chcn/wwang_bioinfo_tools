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
    usage = 'usage: %prog <-g genome_file> <-b bed_file> <-o output_name>'
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
    return optparser


def opt_validate(optparser):
    
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    
    (options, args) = optparser.parse_args()

    # input methylation info file must be given
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

    # output bed file must be given
    if not options.output:
        sys.stdout.write('Error: Output file must be given!\nExit\n')
        optparser.print_help()
        sys.exit(1)
        if os.path.isfile(os.path.expanduser(options.output)):
            sys.stdout.write(f'Warning: Output file: {options.output} alread exist!\n')
            optparser.print_help()
            sys.exit(1)

    return options


def cal_sequence_feature(sequence):
	'''
	Input: sequence
	Output: GC content, CpG density, CpG ratio
	'''
	n = len(sequence)
	G = sequence.upper().count('G')
	C = sequence.upper().count('C')
	CpG = sequence.upper().count('CG')
	return 1.0 * (G + C) / n, 1.0 * CpG / n, 1.0 * CpG * n / C / G


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
         open(os.path.expanduser(options.output)) as output_fhd:
        output_csv = csv.writer(output_fhd)
        output_csv.writerow(['name','GC content','CpG density', 'CpG ratio'])
        for line in input_fhd:
        	if not line:
        		continue
        	line = line.strip().split()
        	chrom, start, end, name = line[:3]
        	start, end = int(start), int(end)
        	start = max(start, 0)
        	sequence = genome[chrom][start:end]
        	GC, CpG_density, CpG_ratio = cal_sequence_feature(sequence)
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