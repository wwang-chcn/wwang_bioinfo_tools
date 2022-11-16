#! /usr/bin/env python3

from collections import defaultdict
from itertools import product


def reverse_complement(seq, complementary_dict):
    return ''.join([complementary_dict[c] for c in seq])


def count_kmer(seq, kmer_len, seq_tpye='DNA'):
    """Count kmer occurance times in sequence"""
    # constant
    alphabet_dict = {'DNA': 'ATCG'}
    alphabet_complementary_dict = {'DNA': {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}}
    # pre-processing
    seq = seq.upper()
    kmer_count_plus = defaultdict(int)
    kmer_count = defaultdict(int)
    # count plus strand
    for i in range(len(seq)-kmer_len):
        kmer_count_plus[seq[i:i+kmer_len]] += 1
    # count minus strand
    for kmer in product(alphabet_dict[seq_tpye], repeat=kmer_len):
        kmer = ''.join(kmer)
        kmer_count[kmer] = kmer_count_plus[kmer] + kmer_count_plus[reverse_complement(kmer, alphabet_complementary_dict[seq_tpye])]
    return kmer_count
