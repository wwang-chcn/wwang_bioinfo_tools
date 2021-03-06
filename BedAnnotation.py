#!/usr/bin/env python3

import os
import sys
from optparse import OptionParser


# ------------------------------------
# Classes
# ------------------------------------

class Bed(dict):

    """Store region"""
    
    def __init__(self):
        pass

    def add_one_region(self, chrom, start, end, name):

        """Add one region to Bed"""

        if not chrom in self:
            self[chrom] = []
        self[chrom].append((start, end, name))
        

# ------------------------------------
# Misc Functions
# ------------------------------------

def match_comparsion(match, start, end, region):
    
    """select cloest element"""
    distance = abs((start + end) / 2 - (region[0] + region[1]) / 2)
    if match and match[0] < distance:
        return match
    else:
        return [distance, region[2]]    # distance, name


def region_annotate(bedRegion, elementRegion):
    
    """annotate one region for one element"""

    match = []
    (chrom, start, end) = bedRegion
    if not chrom in elementRegion:
        return match
    for region in elementRegion[chrom]:
        if region[1] <= start:
            continue
        if region[0] < end and region[1] > start:
            match = match_comparsion(match, start, end, region)
            continue
        if region[0] >= end:
            break
    
    return match


# ------------------------------------
# Sub Functions
# ------------------------------------

def prepare_optparser():

    """
    Prepare optparser object. New options will be added in thisfunction first.
    """

    usage = "USAGE: %prog <-b bed> <-a annotation> [-f annotation_format] [-p promoter] [-e enhancer] [--name2]"
    description = "BedAnnotation -- Annotate a bed file."
    
    # option processor
    optparser = OptionParser(version="%prog 0.1", description=description, usage=usage, add_help_option=False)
    
    # basic setting
    optparser.add_option("-h","--help", action="help", help="Show this help message and exit.")

    # input options
    optparser.add_option("-b", "--bed", dest="bed", type="string",\
                         help="Bed file to be annotated.")
    optparser.add_option("-a", "--ann", dest="ann", type="string",\
                         help="Annotation file.")
    optparser.add_option("-f", "--format", dest="format", type="choice", action='store', choices=['refGene', 'refSeq', 'genePredExt'],\
                         help="Annotation file format. Can be genePredExt, refGene or refSeq.\nDefault is genePredExt.", default="genePredExt")

    # annotate options
    optparser.add_option("-p", "--promoter", dest="promoter", type="int",\
                         help="Promoter range from TSS. Default is 3000.", default=3000)
    optparser.add_option("-e", "--enhancer", dest="enhancer", type="int",\
                         help="Enhancer range from TSS. 0 means not annotate enhancer. Must be larger than promoter if not set to zero. Default is 0.", default=0)
    optparser.add_option("--name2", dest="name2", action="store_true",\
                         help="Whether use name2 for annotation. Only support refGene format now. Default is flase.")

    return optparser

def opt_validate(optparser):

    """Validate options from a OptParser object.

    Ret: Validated options object.
    """

    (options, args) = optparser.parse_args()

    # input bed and anntation file must be given
    if not (options.bed and options.ann):
        sys.stdout.write("Input bed and anntation file must be given!\n")
        optparser.print_help()
        sys.exit(1)

    if not os.path.isfile(options.bed):
        sys.stdout.write("Can not find bed file: {}.\n".format(options.bed))
        optparser.print_help()
        sys.exit(1)

    if not os.path.isfile(options.ann):
        sys.stdout.write("Can not find anntation file: {}.\n".format(options.ann))
        optparser.print_help()
        sys.exit(1)

    # validate promoter and enhancer range
    if options.enhancer != 0:
        if options.enhancer <= options.promoter:
            sys.stdout.write("Enhancer must be larger than promoter if not set to zero.\n")
            optparser.print_help()
            sys.exit(1)

    return options

def load_annotation(anntation, format, promoter, enhancer, name2):

    """Load anntation file & return annotationRegion, TSS"""
    promoterRegion, exonRegion, intronRegion, enhancerRegion = Bed(), Bed(), Bed(), Bed()
    annotationRegion = {"promoter": promoterRegion, "exon": exonRegion, "intron": intronRegion, "enhancer": enhancerRegion}
    TSS = {}

    annFhd = open(anntation)

    if format == 'genePredExt':
        chrom_index = 2
        strand_index = 3
        name_index = 12 if name2 else 1
        for line in annFhd:
            line = line.strip().split()

            # load parameters & TSS
            chrom, strand, name = line[chrom_index-1], line[strand_index-1], line[name_index-1]
            exons = list(zip(line[8].split(","), line[9].split(",")))[:-1]
            introns = list(zip(line[9].split(","), line[8].split(",")[1:]))[:-1]
            if strand == "+":
                tss = int(line[3])
            else:
                tss = int(line[4])
            TSS[name] = (chrom, tss)

            # get region
            promoterRegion.add_one_region(chrom, max(0,tss-promoter), tss+promoter, name)
            for exon in exons:
                exonRegion.add_one_region(chrom, int(exon[0]), int(exon[1]), name)
            for intron in introns:
                intronRegion.add_one_region(chrom, int(intron[0]), int(intron[1]), name)
            enhancerRegion.add_one_region(chrom, tss-enhancer, tss+enhancer, name)

    elif format == "refGene":
        chrom_index = 3
        strand_index = 4
        name_index = 13 if name2 else 2
        for line in annFhd:
            line = line.strip().split()

            # load parameters & TSS
            chrom, strand, name = line[chrom_index-1], line[strand_index-1], line[name_index-1]
            exons = list(zip(line[9].split(","), line[10].split(",")))[:-1]
            introns = list(zip(line[10].split(","), line[9].split(",")[1:]))[:-1]
            if strand == "+":
                tss = int(line[4])
            else:
                tss = int(line[5])
            TSS[name] = (chrom, tss)

            # get region
            promoterRegion.add_one_region(chrom, max(0,tss-promoter), tss+promoter, name)
            for exon in exons:
                exonRegion.add_one_region(chrom, int(exon[0]), int(exon[1]), name)
            for intron in introns:
                intronRegion.add_one_region(chrom, int(intron[0]), int(intron[1]), name)
            enhancerRegion.add_one_region(chrom, tss-enhancer, tss+enhancer, name)

    elif format == "refSeq":
        for line in annFhd:
            if not line:
                continue
            line = line.strip().split()

            # load parameters & TSS
            chrom, strand, name, txStart, exonNum = line[0], line[5], line[3], int(line[1]), int(line[9])
            lengths = [int(x) for x in line[10].split(',')]
            starts = [int(x) for x in line[11].split(',')]
            if strand == '+':
                tss = int(line[1])
            else:
                tss = int(line[2])
            TSS[name] = (chrom, tss)

            # get region
            promoterRegion.add_one_region(chrom, max(0,tss-promoter), tss+promoter, name)
            try:    # intron number is 1 less than exon number
                for i in range(exonNum):
                    exonRegion.add_one_region(chrom, txStart+starts[i], txStart+starts[i]+lengths[i], name)
                    intronRegion.add_one_region(chrom, txStart+starts[i]+lengths[i], txStart+starts[i+1], name)
            except IndexError:
                pass
            enhancerRegion.add_one_region(chrom, tss-enhancer, tss+enhancer, name)

    #for region in annotationRegion.itervalues():
    #    for chrom in region.keys():
    #        region[chrom].sort(key = lambda x: x[0])
    for label in annotationRegion:
        for chrom in annotationRegion[label]:
            annotationRegion[label][chrom].sort(key = lambda x: x[0])
    annFhd.close()

    return (annotationRegion, TSS)


# ------------------------------------
# Main function
# ------------------------------------

def main():

    # read the options and validate them
    options = opt_validate(prepare_optparser())

    # load annotation information & TSS
    (annotationRegion, TSS) = load_annotation(options.ann, options.format, options.promoter, options.enhancer, options.name2)

    # do annotation
    bedFhd = open(options.bed)
    rFhd = open(options.bed.rsplit(".",1)[0]+"_annotation.txt", "w")
    for line in bedFhd:
        rFhd.write(line.strip()+"\t")
        line = line.strip().split()
        chrom, start, end = line[0], int(line[1]), int(line[2])
        if not chrom in annotationRegion['promoter']:
            rFhd.write("intergenic\tNA\tNA\n")
            continue
        for element in ["promoter", "exon", "intron", "enhancer"]:    # priority
            match = region_annotate((chrom, start, end), annotationRegion[element])
            if any(match):
                rFhd.write("{}\t{}\t{}\n".format(element, match[1], int(match[0])))
                break
        else:
            rFhd.write("intergenic\tNA\tNA\n")
    bedFhd.close()
    rFhd.close()

# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write("User interrupts me! ;-) See you ^.^!")
        sys.exit(0)
