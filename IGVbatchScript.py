#!/usr/bin/env python

# migrate to python3
# Apr-13-2016

import os, sys
from optparse import OptionParser, OptionGroup


# ------------------------------------
# Sub Functions
# ------------------------------------
def prepare_optparser():
    
    """
    Prepare optparser object. New options will be added in thisfunction first.
    """
    
    program_name = os.path.basename(sys.argv[0])
    usage = f'USAGE: {program_name} <-n name> <-b bed> <-d dir> [-m method]'
    description = 'IGVbatchScript -- generate IGV batch Script.'
    
    # option processor
    optparser = OptionParser(version=f'{program_name} 0.1', description=description, usage=usage, add_help_option=False)
    
    # basic setting
    optparser.add_option('-h', '--help', action='help', help='Show this help message and exit.')
    optparser.add_option('-n', '--name', dest='name', type='string',
                         help='Name of this run. If not given, the time will be used. The output files will in the directory of this name under the running directory.')
    optparser.add_option('-b', '--bed', dest='bed', type='string',\
                         help='Bed region to be capture.')
    optparser.add_option('-d', '--dir', dest='dir', type='string',\
                         help='Directory to be store snapshots.')
    optparser.add_option('-m', '--method', dest='method', type='string',\
                         help='Method use for naming snapshots. Can be "ordinal" or "region". Default is ordinal.', default='ordinal')
    
    return optparser

def opt_validate(optparser):
    
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # input name, bed, dir must be given
    if not (options.name and options.bed and options.dir):
        print('input name, bed, dir must be given!\n')
        optparser.print_help()
        sys.exit(1)
    
    # input bed and dir must valid
    if not os.path.isfile(options.bed):
        print('illegal bed options!\n')
        optparser.print_help()
        sys.exit(1)
    
    if options.dir == '/':
        options.dir = parentPath[:-1]
    parentPath = options.dir
    parentPath = parentPath.rsplit('/',1)[0]
    if not os.path.isdir(parentPath):
        sys.stdout(f'Warning! Directory {parentPath} is not exist.\n')
    
    # input method must be valid
    if options.method:
        if not (options.method == 'ordinal' or options.method == 'region'):
            print('method is not valid')
            optparser.print_help()
            sys.exit(1)
    
    return options


# ------------------------------------
# Main function
# ------------------------------------
def main():
    
    # read the options and validate them
    options = opt_validate(prepare_optparser())
    
    if not os.path.isdir(options.dir):
        os.mkdir(os.path.abspath(options.dir))
    
    rfhd = open(options.name+'.txt', 'w')
    rfhd.write(f'snapshotDirectory {os.path.abspath(options.dir)}\n')
    bedfhd = open(options.bed)
    i = 0
    for line in bedfhd:
        if not line:
            continue
        i += 1
        line = line.strip().split('\t')
        if len(line) < 3:
            print(f'Line {i} in {options.bed} has less than 3 fileds!\n')
            sys.exit(1)
        if options.method == 'ordinal':
            rfhd.write(f'goto {line[0]}:{line[1]}-{line[2]}\nsnapshot {options.name}_{i}.png\n')
        if options.method == 'region':
            rfhd.write(f'goto {line[0]}:{line[1]}-{line[2]}\nsnapshot {options.name}_{line[0]}_{line[1]}_{line[2]}.png\n')
    rfhd.close()
    bedfhd.close()


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn('User interrupts me! ;-) See you ^.^!')
        sys.exit(0)