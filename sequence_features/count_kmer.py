#! /usr/bin/env python3

from collections import defaultdict
from itertools import product

# constant
alphabet_dict = {'DNA': 'ATCGN'}
alphabet_complementary_dict = {'DNA': {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}}


def reverse_complement(seq, complementary_dict):
    return ''.join([complementary_dict[c] for c in seq[::-1]])


def count_kmer(seq, kmer_len, seq_tpye='DNA', reverse_complement=True):
    """Count k-mer occurance times in sequence."""
    # pre-processing
    seq = seq.upper()
    kmer_count_plus = defaultdict(int)
    kmer_count = defaultdict(int)
    # count plus strand
    for i in range(len(seq)-kmer_len+1):
        kmer_count_plus[seq[i:i+kmer_len]] += 1
    if reverse_complement: # count minus strand
        for kmer in product(alphabet_dict[seq_tpye], repeat=kmer_len):
            kmer = ''.join(kmer)
            kmer_count[kmer] = kmer_count_plus[kmer] + kmer_count_plus[reverse_complement(kmer, alphabet_complementary_dict[seq_tpye])]
        return kmer_count
    else:
        return kmer_count_plus


def count_kmer2(seq, kmer_len, seq_tpye='DNA', reverse_complement=True):
    """Count the occurance times k-mer from length 1 to k in sequence."""
    # pre-processing
    seq = seq.upper()
    kmer_count_plus = defaultdict(int)
    kmer_count = defaultdict(int)
    # count plus strand
    #     initial stage
    for i in range(min(len(seq), kmer_len)):
        for j in range(1,i+1):
            kmer_count_plus[seq[i-j:i]] += 1
        print(i, kmer_count_plus)
    #     whole sequence
    for i in range(len(seq)-kmer_len+1):
        for j in range(kmer_len):
            kmer_count_plus[seq[i+j:i+kmer_len]] += 1
        print(i+kmer_len, kmer_count_plus)
    if reverse_complement: # count minus strand
        for i in range(1, kmer_len + 1):
            for kmer in product(alphabet_dict[seq_tpye], repeat=i):
                kmer = ''.join(kmer)
                kmer_count[kmer] = kmer_count_plus[kmer] + kmer_count_plus[reverse_complement(kmer, alphabet_complementary_dict[seq_tpye])]
        return kmer_count
    else:
        return kmer_count_plus
