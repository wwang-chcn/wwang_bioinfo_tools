#! /usr/bin/env python3

from collections import defaultdict
from itertools import product

# constant
alphabet_dict = {'DNA': 'ATCGN'}
alphabet_complementary_dict = {'DNA': {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}}


def reverse_complement(seq, complementary_dict):
    return ''.join([complementary_dict[c] for c in seq[::-1]])


def count_kmer(seq, kmer_len, seq_tpye='DNA'):
    """Count k-mer occurance times in sequence."""
    # pre-processing
    seq = seq.upper()
    kmer_count_plus = defaultdict(int)
    kmer_count = defaultdict(int)
    # count plus strand
    for i in range(len(seq)-kmer_len+1):
        kmer_count_plus[seq[i:i+kmer_len]] += 1
    # count minus strand
    for kmer in product(alphabet_dict[seq_tpye], repeat=kmer_len):
        kmer = ''.join(kmer)
        kmer_count[kmer] = kmer_count_plus[kmer] + kmer_count_plus[reverse_complement(kmer, alphabet_complementary_dict[seq_tpye])]
    return kmer_count

