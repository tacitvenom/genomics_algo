import math
import numpy as np
from collections import Counter
from typing import Sequence

from numpy.core.multiarray import concatenate

COMPLEMENTARY_BASE = {"A": "T", "C": "G", "G": "C", "T": "A"}


def longest_common_prefix(s1, s2):
    """
    Finds the longest common prefix (substring) given two strings

    s1: string
        First string to compare
    s2: string
        Second string to compare

    Returns:
    longest_common_prefix: String
        Longest common prefix between s1 and s2
    """
    i = 0
    while i < min(len(s1), len(s2)):
        if s1[i] != s2[i]:
            break
        i += 1
    return s1[:i]


def reverse_complement(s):
    """
    Find the reverse complement of a DNA strand
    s: string
        A DNA sequence of a strand - the string must have the characters A, C, G and T

    Returns:
        DNA sequence of the opposite strand in the reverse order
    """
    # TODO: make robust againt garbage values
    return "".join([COMPLEMENTARY_BASE[base] for base in s[::-1]])


def read_genome(filename):
    """
    Reads a genome from a .fa file

    filename: string
        relative or absolute path of the .fa file to be read from

    Returns:
    genome: string
        Genome string
    """
    with open(filename) as f:
        genome = "".join(
            [line.rstrip() for line in f.readlines() if not line.startswith(">")]
        )
    return genome


def read_fastq(filename):
    """
    Reads sequences and qualities from a .fastq file

    filename: string
        relative or absolute path of the .fa file to be read from

    Returns:
    reads: list
        List of sequence reads
    qualities: list
        List of qualities corresponding to each sequence read

    """
    reads = []
    qualities = []
    with open(filename, "r") as f:
        while True:
            f.readline()
            read = f.readline().rstrip()
            f.readline()
            seq_qualities = f.readline().rstrip()
            if len(read) == 0:
                break
            reads.append(read)
            qualities.append(seq_qualities)

    return reads, qualities


def map_phred33_to_error_probability(phred33):
    """Maps a ASCII phred33 quality character to a quality score in fraction
    # TODO
    Args:
        quality ([type]): [description]
    """
    return 10 ** (-phred33 / 10)


def map_errorprobability_to_phred33(probability):
    """
    # TODO
    """
    return -10 * math.log(probability, 10)


def map_phred33_ascii_to_qualityscore(phred33_char):
    """Maps a ASCII phred33 quality character to a quality score
    # TODO
    Args:
        quality ([type]): [description]
    """
    return ord(phred33_char) - 33


def get_freq_for_qualities(qualities):
    """
    Generates a frequency distribution from a list of quality strings
    """
    concatenated_qualities = ''.join(qualities)
    quality_scores = [map_phred33_ascii_to_qualityscore(char) for char in concatenated_qualities] 
    freq = Counter(quality_scores)
    sorted_freq = sorted(dict(freq).items())
    return [pair[0] for pair in sorted_freq], [pair[1] for pair in sorted_freq]


def same_length_reads(reads):
    """
    Returns true if the list of sequencing reads has at least one sequence read and
    each of the reads has the same length, False otherwise
    """
    assert len(reads) > 0
    lengths = [len(read) for read in reads]
    return min(lengths) == max(lengths)


def find_GC_by_position(reads):
    """
    Returns the average GC content per index in a list of sequencing reads
    """
    assert same_length_reads(reads)
    reads_length = len(reads[0])
    gc = np.zeros(reads_length)
    for read in reads:
        for index, base in enumerate(read):
            if base in ['G', 'C']:
                gc[index] += 1
    gc = gc/len(reads)
    return gc


def get_base_freq(reads):
    """
    Returns the aggregate frequency of bases in the sequencing reads
    """
    concatenated_reads = ''.join(reads)
    base_freq = Counter(concatenated_reads)
    return base_freq
