#!/usr/bin/env python3

import os
import sys
from collections import defaultdict
from optparse import OptionParser
from typing import Optional


# ------------------------------------
# Classes
# ------------------------------------
class Bed():
    """Store region"""

    def __init__(self):
        self.regions = defaultdict(list)

    def add_one_region(self, chrom: str, start: int, end: int,
                       name: str) -> None:
        """Add one region to Bed"""

        self.regions[chrom].append((start, end, name))

    def __contains__(self, chrom: str) -> bool:
        return self.regions.__contains__(chrom)

    def items(self):
        return self.regions.items()

    def __getitem__(self, key):
        return self.regions.__getitem__(key)


# ------------------------------------
# Misc Functions
# ------------------------------------
def match_comparsion(match: Optional[list], start: int, end: int, region):
    """Compare new match with previous one, and select the most cloest element

Parameters
----------
match : list[distance, name] or None
start : int
end : int
region : list[start: int, end: int, ...]

Returns
----------
match : list[distance, name]
status: bool
        If the new match is far than previous one.
        If true, the target annotation have been found.
"""
    distance = abs((start + end) - (region[0] + region[1])) >> 1
    if match is None:
        return [distance, region[2]], False
    if match[0] < distance:
        return match, True
    else:
        return [distance, region[2]], False


def region_annotate(bedRegion, elementRegion):
    """annotate one region for one element"""

    match = None
    (chrom, start, end) = bedRegion
    if chrom not in elementRegion:
        return match
    for region in elementRegion[chrom]:
        if region[1] <= start:
            continue
        if region[0] < end and region[1] > start:
            match, status = match_comparsion(match, start, end, region)
            if status:
                break
            else:
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
    optparser = OptionParser(version="%prog 0.1",
                             description=description,
                             usage=usage,
                             add_help_option=False)

    # basic setting
    optparser.add_option("-h",
                         "--help",
                         action="help",
                         help="Show this help message and exit.")

    # input options
    optparser.add_option("-b",
                         "--bed",
                         dest="bed",
                         type="string",
                         help="Bed file to be annotated.")
    optparser.add_option("-a",
                         "--ann",
                         dest="ann",
                         type="string",
                         help="Annotation file.")
    optparser.add_option(
        "-f",
        "--format",
        dest="format",
        type="choice",
        action='store',
        choices=['refGene', 'refSeq', 'genePredExt'],
        help="Annotation file format. Can be genePredExt, refGene or refSeq.\nDefault is genePredExt.",
        default="genePredExt")

    # annotate options
    optparser.add_option("-p",
                         "--promoter",
                         dest="promoter",
                         type="int",
                         help="Promoter range from TSS. Default is 3000.",
                         default=3000)
    optparser.add_option(
        "-e",
        "--enhancer",
        dest="enhancer",
        type="int",
        help="Enhancer range from TSS. 0 means not annotate enhancer. Must be larger than promoter if not set to zero. Default is 0.",
        default=0)
    optparser.add_option(
        "--name2",
        dest="name2",
        action="store_true",
        help="Whether use name2 for annotation. Only support refGene format now. Default is flase."
    )

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
        sys.stdout.write("Can not find anntation file: {}.\n".format(
            options.ann))
        optparser.print_help()
        sys.exit(1)

    # validate promoter and enhancer range
    if options.enhancer != 0:
        if options.enhancer <= options.promoter:
            sys.stdout.write(
                "Enhancer must be larger than promoter if not set to zero.\n")
            optparser.print_help()
            sys.exit(1)

    return options


def load_annotation(anntation, format, promoter, enhancer, name2):
    """Load anntation file & return annotationRegion, TSS"""
    promoterRegion, exonRegion, intronRegion, enhancerRegion = Bed(), Bed(
    ), Bed(), Bed()
    annotationRegion = {
        "promoter": promoterRegion,
        "exon": exonRegion,
        "intron": intronRegion,
        "enhancer": enhancerRegion
    }
    TSS = {}

    annFhd = open(anntation)

    if format == 'genePredExt':
        chrom_index = 2
        strand_index = 3
        name_index = 12 if name2 else 1
        for line in annFhd:
            line = line.strip().split()

            # load parameters & TSS
            chrom, strand, name = line[chrom_index -
                                       1], line[strand_index -
                                                1], line[name_index - 1]
            exons = list(zip(line[8].split(","), line[9].split(",")))[:-1]
            introns = list(zip(line[9].split(","),
                               line[8].split(",")[1:]))[:-1]
            if strand == "+":
                tss = int(line[3])
            else:
                tss = int(line[4])
            TSS[name] = (chrom, tss)

            # get region
            promoterRegion.add_one_region(chrom, max(0, tss - promoter),
                                          tss + promoter, name)
            for exon in exons:
                exonRegion.add_one_region(chrom, int(exon[0]), int(exon[1]),
                                          name)
            for intron in introns:
                intronRegion.add_one_region(chrom, int(intron[0]),
                                            int(intron[1]), name)
            enhancerRegion.add_one_region(chrom, tss - enhancer,
                                          tss + enhancer, name)

    elif format == "refGene":
        chrom_index = 3
        strand_index = 4
        name_index = 13 if name2 else 2
        for line in annFhd:
            line = line.strip().split()

            # load parameters & TSS
            chrom, strand, name = line[chrom_index -
                                       1], line[strand_index -
                                                1], line[name_index - 1]
            exons = list(zip(line[9].split(","), line[10].split(",")))[:-1]
            introns = list(zip(line[10].split(","),
                               line[9].split(",")[1:]))[:-1]
            if strand == "+":
                tss = int(line[4])
            else:
                tss = int(line[5])
            TSS[name] = (chrom, tss)

            # get region
            promoterRegion.add_one_region(chrom, max(0, tss - promoter),
                                          tss + promoter, name)
            for exon in exons:
                exonRegion.add_one_region(chrom, int(exon[0]), int(exon[1]),
                                          name)
            for intron in introns:
                intronRegion.add_one_region(chrom, int(intron[0]),
                                            int(intron[1]), name)
            enhancerRegion.add_one_region(chrom, tss - enhancer,
                                          tss + enhancer, name)

    elif format == "refSeq":
        for line in annFhd:
            if not line:
                continue
            line = line.strip().split()

            # load parameters & TSS
            chrom, strand, name, txStart, exonNum = line[0], line[5], line[
                3], int(line[1]), int(line[9])
            lengths = [int(x) for x in line[10].split(',')]
            starts = [int(x) for x in line[11].split(',')]
            if strand == '+':
                tss = int(line[1])
            else:
                tss = int(line[2])
            TSS[name] = (chrom, tss)

            # get region
            promoterRegion.add_one_region(chrom, max(0, tss - promoter),
                                          tss + promoter, name)
            try:  # intron number is 1 less than exon number
                for i in range(exonNum):
                    exonRegion.add_one_region(chrom, txStart + starts[i],
                                              txStart + starts[i] + lengths[i],
                                              name)
                    intronRegion.add_one_region(
                        chrom, txStart + starts[i] + lengths[i],
                        txStart + starts[i + 1], name)
            except IndexError:
                pass
            enhancerRegion.add_one_region(chrom, tss - enhancer,
                                          tss + enhancer, name)

    for label in annotationRegion:
        for chrom, values in annotationRegion[label].items():
            values.sort(key=lambda x: x[0])
    annFhd.close()

    return (annotationRegion, TSS)


# ------------------------------------
# Main function
# ------------------------------------
def main():

    # read the options and validate them
    options = opt_validate(prepare_optparser())

    # load annotation information & TSS
    (annotationRegion, _) = load_annotation(options.ann, options.format,
                                            options.promoter, options.enhancer,
                                            options.name2)

    # do annotation
    output_file_name = f'{options.bed.rsplit(".", 1)[0]}_annotation.txt'
    no_match_info = '\tintergenic\tNA\tNA\n'
    elsments_priority = ['promoter', 'exon', 'intron', 'enhancer']

    with open(options.bed) as bed_fhd, \
         open(output_file_name, 'w') as output_fhd:
        for line in bed_fhd:
            new_line = line.strip()
            chrom, start, end, *_ = line.strip().split()
            start, end = int(start), int(end)
            if chrom not in annotationRegion['promoter']:
                new_line += no_match_info
                continue
            for element in elsments_priority:
                match = region_annotate((chrom, start, end),
                                        annotationRegion[element])
                if match is not None:
                    new_line += f'\t{element}\t{match[1]:s}\t{match[0]:d}\n'
                    break
            else:
                new_line += no_match_info
            output_fhd.write(new_line)


# ------------------------------------
# Program running
# ------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stdout.write("User interrupts me! ;-) See you ^.^!")
        sys.exit(0)
