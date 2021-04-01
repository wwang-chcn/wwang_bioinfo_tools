#! /usr/bin/env python3

# Jan-18-2019

import os, sys
import pysam
from collections import defaultdict

def alignment_reads_snp_count(line,snp,quality_score_threshold=30):
    strain1, strain2, match, mismatch = 0, 0, 0, 0
    for align in line.get_aligned_pairs(matches_only=True):
        if not align[1] in snp[line.reference_name]: # no snp in this position
            continue
        if ord(line.qual[align[0]])-33 < quality_score_threshold:
            continue
        if line.seq[align[0]] == snp[line.reference_name][align[1]][0]:
            strain1 += 1
            match += 1
        elif line.seq[align[0]] == snp[line.reference_name][align[1]][1]:
            strain2 += 1
            match += 1
        else:
            mismatch += 1
    return strain1, strain2, match, mismatch


def bam_split(filterSNP,inputBAM,outputBAM1,outputBAM2,quality_score_threshold=30):
    def calculate_plines(pre_lines, outputFhd1, outputFhd2, snp, statistics, quality_score_threshold=30):
        '''statistics = [match, mismatch, strain1, strain2, s1_reads, s2_reads, mixed, nosnp]
        '''
        match, mismatch, strain1, strain2 = 0, 0, 0, 0
        for line in pre_lines:
            if line.is_unmapped or not line.reference_name.startswith("chr") or '_' in line.reference_name:
                continue
            if not line.reference_name in snp:
                continue
            strain1_, strain2_, match_, mismatch_ = alignment_reads_snp_count(line,snp,quality_score_threshold)
            match += match_
            mismatch += mismatch_
            strain1 += strain1_
            strain2 += strain2_
        statistics[0] += match
        statistics[1] += mismatch
        statistics[2] += strain1
        statistics[3] += strain2
        if strain1 + strain2 < 2: # no snp
            statistics[7] += 1
        elif 1.0 * strain1 / (strain1 + strain2) >= 2.0 / 3: # s1 reads group
            statistics[4] += 1
            for line in pre_lines:
                outputFhd1.write(line)
        elif 1.0 * strain2 / (strain1 + strain2) >= 2.0 / 3: # s2 reads group
            statistics[5] += 1
            for line in pre_lines:
                outputFhd2.write(line)
        else:
            statistics[6] += 1
    # load snp
    snp = defaultdict(dict)
    with open(filterSNP) as fhd:
        for line in fhd:
            chrom, start, end, n1, n2, *_ = line.strip().split()
            snp[chrom][int(start)] = [n1,n2]
    total = 0
    statistics = [0] * 8 # match, mismatch, strain1, strain2, s1_reads, s2_reads, mixed, nosnp
    with pysam.Samfile(inputBAM,'rb') as inputFhd, \
         pysam.Samfile(outputBAM1,'wb',template=inputFhd) as outputFhd1, \
         pysam.Samfile(outputBAM2,'wb',template=inputFhd) as outputFhd2:
        pre_lines = []
        for line in inputFhd:
            if pre_lines and pre_lines[0].query_name == line.query_name: # same reads group
                pre_lines.append(line)
            else: # new reads group
                if pre_lines:
                    total += 1
                    calculate_plines(pre_lines, outputFhd1, outputFhd2, snp, statistics,quality_score_threshold)
                pre_lines = [line]
        if pre_lines:
            total += 1
            calculate_plines(pre_lines, outputFhd1, outputFhd2, snp, statistics,quality_score_threshold)
    statistics.append(total)
    return statistics

def alignment_reads_snp_count_inspect(line,snp,quality_score_threshold=30):
    strain1, strain2, match, mismatch = 0, 0, 0, 0
    for align in line.get_aligned_pairs(matches_only=True,with_seq=True):
        if not align[1] in snp[line.reference_name]: # no snp in this position
            if align[2].upper() != line.seq[align[0]].upper():
                sys.stdout.write(f'chr: {line.reference_name} pos: {align[1]}. No snp in this position. ref: {align[2]} query: {line.seq[align[0]]}\n')
            continue
        if ord(line.qual[align[0]])-33 < quality_score_threshold:
            sys.stdout.write(f'Low sequencing position quality. ref: {snp[line.reference_name][align[1]][0]}, {snp[line.reference_name][align[1]][1]} query: {line.seq[align[0]]}\n')
            continue
        if line.seq[align[0]] == snp[line.reference_name][align[1]][0]:
            sys.stdout.write(f'chr: {line.reference_name} pos: {align[1]}. Match. ref: {snp[line.reference_name][align[1]][0]}, {snp[line.reference_name][align[1]][1]} query: {line.seq[align[0]]}\n')
            strain1 += 1
            match += 1
        elif line.seq[align[0]] == snp[line.reference_name][align[1]][1]:
            sys.stdout.write(f'chr: {line.reference_name} pos: {align[1]}. Match. ref: {snp[line.reference_name][align[1]][0]}, {snp[line.reference_name][align[1]][1]} query: {line.seq[align[0]]}\n')
            strain2 += 1
            match += 1
        else:
            sys.stdout.write(f'chr: {line.reference_name} pos: {align[1]}. Misatch. ref: {snp[line.reference_name][align[1]][0]}, {snp[line.reference_name][align[1]][1]} query: {line.seq[align[0]]}\n')
            mismatch += 1
    return strain1, strain2, match, mismatch

def bam_inspect(filterSNP,inputBAM,quality_score_threshold=30):
    def calculate_plines(pre_lines, snp, quality_score_threshold=30):
        match, mismatch, strain1, strain2 = 0, 0, 0, 0
        for line in pre_lines:
            sys.stdout.write(f'{line.to_string()}\n')
            if line.is_unmapped or not line.reference_name.startswith("chr") or '_' in line.reference_name:
                sys.stdout.write('Umapped or small chromsome\n')
                continue
            if not line.reference_name in snp:
                sys.stdout.write('Non-snp chromsome\n')
                continue
            strain1_, strain2_, match_, mismatch_ = alignment_reads_snp_count_inspect(line,snp,quality_score_threshold)
            match += match_
            mismatch += mismatch_
            strain1 += strain1_
            strain2 += strain2_
        if strain1 + strain2 < 2: # no snp
            sys.stdout.write('No SNP\n')
        elif 1.0 * strain1 / (strain1 + strain2) >= 2.0 / 3: # s1 reads group
            sys.stdout.write('Strain1\n')
        elif 1.0 * strain2 / (strain1 + strain2) >= 2.0 / 3: # s2 reads group
            sys.stdout.write('Strain2\n')
        else:
            sys.stdout.write('Mixed\n')
    # load snp
    snp = defaultdict(dict)
    with open(filterSNP) as fhd:
        for line in fhd:
            chrom, start, end, n1, n2, *_ = line.strip().split()
            snp[chrom][int(start)] = [n1,n2]
    with pysam.Samfile(inputBAM,'rb') as inputFhd:
        pre_lines = []
        for line in inputFhd:
            if pre_lines and pre_lines[0].query_name == line.query_name: # same reads group
                pre_lines.append(line)
            else: # new reads group
                if pre_lines:
                    calculate_plines(pre_lines, snp)
                    yield
                pre_lines = [line]
        if pre_lines:
            calculate_plines(pre_lines, snp)
            yield

# ----- parameters -----
filterSNP = sys.argv[1]
inputBAM = sys.argv[2]
outputBAM1 = sys.argv[3]
outputBAM2 = sys.argv[4]

match, mismatch, strain1, strain2, s1_reads, s2_reads, mixed, nosnp, total = bam_split(filterSNP,inputBAM,outputBAM1,outputBAM2)
output_template = f'''\
Total: {total:,}
    match: {match:,}
        strain1: {strain1:,}
        strain2: {strain2:,}
    mismatch: {mismatch:,}

    strain1 reads group: {s1_reads:,}
    strain2 reads group: {s2_reads:,}
    mixed: {mixed:,}
    nosnp: {nosnp:,}
'''
sys.stdout.write(output_template)
