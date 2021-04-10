#! /usr/bin/env python3

import os, sys
from pysam import AlignmentFile
from optparse import OptionParser
from contextlib import contextmanager


# ------------------------------------
# Sub Functions
# ------------------------------------

def prepare_optparser():
    
    """
    Prepare optparser object. New options will be added in thisfunction first.
    """
    
    usage = "USAGE: {} [-a alignmentFile] [-o outputBedFile]".format(os.path.basename(sys.argv[0]))
    description = "bamToBed -- Get reads cover region from sam/bam file (especially for RNA-seq samples)."
    
    # option processor
    optparser = OptionParser(version="{} 0.1".format(os.path.basename(sys.argv[0])), description=description, usage=usage, add_help_option=False)
    
    # basic setting
    optparser.add_option("-h","--help", action="help", help="Show this help message and exit.")
    optparser.add_option("-a",dest="alignmentFile", type="string",\
                         help="AlignmentFile. Sam/Bam format. Default is stdin.")
    optparser.add_option("-o",dest="outputBedFile", type="string",\
                         help="outputBedFile. Default is stdout.")
    
    return optparser

def opt_validate(optparser):
    
    """Validate options from a OptParser object.
    
    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    if not options.alignmentFile:
        options.alignmentFile = '-'
    else:
        if not os.path.isfile(options.alignmentFile):
            sys.stdout.write("Error! Can not find alignmentFile: {}!\nExit!\n".format(options.alignmentFile))
            sys.exit(1)
    
    if not options.outputBedFile:
        options.outputBedFile = '-'
    
    return options

def mergeBed(Bed):
    if len(Bed) < 2:
        return Bed
    Bed = sorted(Bed)
    new_Bed = []
    chrom, start, end = Bed[0]
    for region in Bed[1:]:
        if region[0] != chrom: # different chromsome
            new_Bed.append((chrom,start,end))
            chrom, start, end = region
            continue
        if region[1] > end: # non-overlap
            new_Bed.append((chrom,start,end))
            chrom, start, end = region
        else: # overlap
            end = max(end, region[2])
    new_Bed.append((chrom,start,end))
    return new_Bed

def mergeAndOutpuBed(Bed,fhd,name,strand="."):
    Bed = mergeBed(Bed)
    if Bed:
        for region in Bed:
            fhd.write("{}\t{}\t{}\t{}\t0\t{}\n".format(region[0],region[1],region[2],name,strand))

@contextmanager
def smart_out_open(filename=None,mode='w'):
    if filename and filename != 'sys.stdout' and filename != '-':
        fhd = open(filename, mode)
    else:
        fhd = sys.stdout
    
    try:
        yield fhd
    finally:
        if fhd is not sys.stdout:
            fhd.close()


def bamToBed(alignmentFile,outputBedFile):
    with AlignmentFile(alignmentFile) as inputFhd, \
        smart_out_open(outputBedFile,"w") as outputFhd:
        regions, name, strand = None, None, None
        for alignmentSegement in inputFhd:
            if alignmentSegement.is_unmapped: # unmapped reads
                continue
            if alignmentSegement.qname != name: # reads from new fragments
                if name:
                    mergeAndOutpuBed(regions,outputFhd,name,strand)
                regions, name, strand = [(alignmentSegement.reference_name,x[0],x[1]) for x in alignmentSegement.get_blocks()], alignmentSegement.qname, "-" if alignmentSegement.is_reverse and alignmentSegement.is_read1 or not alignmentSegement.is_reverse and not alignmentSegement.is_read1 else "+"
            else: # reads from same fragments
                regions.extend([(alignmentSegement.reference_name,x[0],x[1]) for x in alignmentSegement.get_blocks()])
        mergeAndOutpuBed(regions,outputFhd,name,strand)


# ------------------------------------
# Main Functions
# ------------------------------------

def main():
    
    # read the options and validate them
    options = opt_validate(prepare_optparser())
    bamToBed(options.alignmentFile,options.outputBedFile)


# ------------------------------------
# Program running
# ------------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write("User interrupts me! ;-) See you ^.^!\n")
        sys.exit(0)

