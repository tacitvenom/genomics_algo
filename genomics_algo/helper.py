import doctest
import math
import numpy as np
import random

from collections import Counter
from typing import Dict, List, Tuple


class Bases:
    A = "A"
    C = "C"
    G = "G"
    T = "T"


COMPLEMENTARY_BASE = {
    Bases.A: Bases.T,
    Bases.C: Bases.G,
    Bases.G: Bases.C,
    Bases.T: Bases.A,
    "N": "N",
}


def longest_common_prefix(s1: str, s2: str) -> str:
    """
    Finds the longest common prefix (substring) given two strings

    s1: First string to compare
    s2: Second string to compare

    Returns:
        Longest common prefix between s1 and s2

    >>> longest_common_prefix("ACTA", "GCCT")
    ''
    >>> longest_common_prefix("ACTA", "ACT")
    'ACT'
    >>> longest_common_prefix("ACT", "ACTA")
    'ACT'
    >>> longest_common_prefix("GATA", "GAAT")
    'GA'
    >>> longest_common_prefix("ATGA", "")
    ''
    >>> longest_common_prefix("", "GCCT")
    ''
    >>> longest_common_prefix("GCCT", "GCCT")
    'GCCT'
    """
    i = 0
    while i < min(len(s1), len(s2)):
        if s1[i] != s2[i]:
            break
        i += 1
    return s1[:i]


def longest_common_suffix(s1: str, s2: str) -> str:
    """
    Finds the longest common suffix (substring) given two strings

    s1: First string to compare
    s2: Second string to compare

    Returns:
        Longest common suffix between s1 and s2

    >>> longest_common_suffix("ACTA", "GCCT")
    ''
    >>> longest_common_suffix("ACTA", "CTA")
    'CTA'
    >>> longest_common_suffix("CTA", "ACTA")
    'CTA'
    >>> longest_common_suffix("GATAT", "GAATAT")
    'ATAT'
    >>> longest_common_suffix("ACTA", "")
    ''
    >>> longest_common_suffix("", "GCCT")
    ''
    >>> longest_common_suffix("GCCT", "GCCT")
    'GCCT'
    """
    return longest_common_prefix(s1[::-1], s2[::-1])[::-1]


def reverse_complement(s: str) -> str:
    """
    Find the reverse complement of a DNA strand
    s: A DNA sequence of a strand - the string must have the characters A, C, G and T

    Returns:
        DNA sequence of the opposite strand in the reverse order

    >>> reverse_complement("ATGC")
    'GCAT'
    >>> reverse_complement("")
    ''
    """
    # TODO: make robust against garbage values
    return "".join([COMPLEMENTARY_BASE[base] for base in s[::-1]])


def read_genome(filename: str) -> str:
    """
    Reads a genome from a .fa file

    filename: relative or absolute path of the .fa file to be read from

    Returns:
        Genome string
    """
    with open(filename) as f:
        genome = "".join(
            [line.rstrip() for line in f.readlines() if not line.startswith(">")]
        )
    return genome


def read_fastq(filename: str) -> Tuple[List[str], List[str]]:
    """
    Reads sequences and qualities from a .fastq file

    filename: relative or absolute path of the .fa file to be read from

    Returns:
        List of sequence reads
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


def map_phred33_to_error_probability(phred33: str) -> float:
    """Maps a ASCII phred33 quality character to error probability"""
    return 10 ** (-phred33 / 10)


def map_errorprobability_to_phred33(probability: float) -> str:
    """Maps an error probability value to ASCII phred33 quality character"""
    return -10 * math.log(probability, 10)


def map_phred33_ascii_to_qualityscore(phred33_char: str) -> float:
    """Maps a ASCII phred33 quality character to a quality score
    >>> map_phred33_ascii_to_qualityscore("#")
    2
    >>> map_phred33_ascii_to_qualityscore("J")
    41
    """
    return ord(phred33_char) - 33


def get_freq_for_qualities(qualities: List[str]) -> Tuple[List[str], List[int]]:
    """Generates a frequency distribution from a list of quality strings"""
    concatenated_qualities = "".join(qualities)
    quality_scores = [
        map_phred33_ascii_to_qualityscore(char) for char in concatenated_qualities
    ]
    freq = Counter(quality_scores)
    sorted_freq = sorted(dict(freq).items())
    values = [pair[0] for pair in sorted_freq]
    frequencies = [pair[1] for pair in sorted_freq]
    return values, frequencies


def same_length_reads(reads: List[str]) -> bool:
    """
    Returns true if the list of sequencing reads has at least one sequence read and
    each of the reads has the same length, False otherwise
    >>> same_length_reads(["AACGTTA"])
    True
    >>> same_length_reads(["AACGTTA", "CGCGTTT"])
    True
    >>> same_length_reads(["AACGTTA", "CGCGTTT", "GTTAC"])
    False
    """
    assert len(reads) > 0
    lengths = [len(read) for read in reads]
    return min(lengths) == max(lengths)


def find_GC_by_position(reads: List[str]) -> np.ndarray:
    """
    Returns the average GC content per index in a list of sequencing reads
    """
    assert same_length_reads(reads)
    reads_length = len(reads[0])
    gc = np.zeros(reads_length)
    for read in reads:
        assert set(read) - {"A", "C", "G", "T", "N"} == set()
        for index, base in enumerate(read):
            if base in ["G", "C"]:
                gc[index] += 1
    gc = gc / len(reads)
    return gc


def get_base_freq(reads: List[str]):
    """
    Returns the aggregate frequency of bases in the sequencing reads
    >>> get_base_freq(["NAACGTTA"])
    Counter({'A': 3, 'T': 2, 'N': 1, 'C': 1, 'G': 1})
    >>> get_base_freq(["AACGTTA", "CGCGTTT"])
    Counter({'T': 5, 'A': 3, 'C': 3, 'G': 3})
    """
    concatenated_reads = "".join(reads)
    return Counter(concatenated_reads)


def generate_artificial_reads(
    genome: str, number_of_reads: int, read_length: int
) -> List[str]:
    """Generate a set of reads randomly from a genome"""
    reads = []
    for _ in range(number_of_reads):
        start_position = random.randint(0, len(genome) - read_length + 1)
        reads.append(genome[start_position : start_position + read_length + 1])
    return reads


def get_frequency_map(text: str, substring_length: int) -> Dict[str, int]:
    """
    Find the frequency of all substring of length in a given text
    >>> get_frequency_map("GTACGTACC", 1)
    {'G': 2, 'T': 2, 'A': 2, 'C': 3}
    >>> get_frequency_map("GTACGTACC", 2)
    {'GT': 2, 'TA': 2, 'AC': 2, 'CG': 1, 'CC': 1}
    >>> get_frequency_map("GTACGTACC", 4)
    {'GTAC': 2, 'TACG': 1, 'ACGT': 1, 'CGTA': 1, 'TACC': 1}
    >>> get_frequency_map("GTACGTACC", 6)
    {'GTACGT': 1, 'TACGTA': 1, 'ACGTAC': 1, 'CGTACC': 1}
    """
    assert substring_length > 0
    assert len(text) > 0

    freq_map = {}
    for index in range(len(text) - substring_length + 1):
        substr = text[index : index + substring_length]
        if substr in freq_map:
            freq_map[substr] += 1
        else:
            freq_map[substr] = 1
    return freq_map


def find_hamming_distance(s1: str, s2: str) -> int:
    """Compute the Hamming distance between two strings of equal length
    >>> find_hamming_distance("ATG", "ATC")
    1
    >>> find_hamming_distance("ATG", "TGA")
    3
    >>> find_hamming_distance("A", "A")
    0
    >>> find_hamming_distance("ATG", "ATG")
    0
    >>> find_hamming_distance("", "")
    0
    >>> find_hamming_distance("GAGGTAGCGGCGTTTAAC", "GTGGTAACGGGGTTTAAC")
    3
    """
    assert len(s1) == len(s2)
    return sum(1 for i in range(len(s1)) if s1[i] != s2[i])


def find_levenshtein_distance(s1: str, s2: str) -> int:
    """Compute the Levenshtein distance between two strings (i.e., minimum number
    of edits including substitution, insertion and deletion needed in a string to
    turn it into another)
    >>> find_levenshtein_distance("AT", "")
    2
    >>> find_levenshtein_distance("AT", "ATC")
    1
    >>> find_levenshtein_distance("ATG", "ATC")
    1
    >>> find_levenshtein_distance("ATG", "TGA")
    2
    >>> find_levenshtein_distance("ATG", "ATG")
    0
    >>> find_levenshtein_distance("", "")
    0
    >>> find_levenshtein_distance("GAGGTAGCGGCGTTTAAC", "GTGGTAACGGGGTTTAAC")
    3
    >>> find_levenshtein_distance("TGGCCGCGCAAAAACAGC", "TGACCGCGCAAAACAGC")
    2
    >>> find_levenshtein_distance("GCGTATGCGGCTAACGC", "GCTATGCGGCTATACGC")
    2
    """
    # initializing a matrix for with `len(s1) + 1` rows and `len(s2) + 1` columns
    D = [[0 for x in range(len(s2) + 1)] for y in range(len(s1) + 1)]

    # fill first column
    for i in range(len(s1) + 1):
        D[i][0] = i
    # fill first row
    for j in range(len(s2) + 1):
        D[0][j] = j

    # fill rest of the matrix
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            distLeft = D[i][j - 1] + 1
            distAbove = D[i - 1][j] + 1
            # import pdb; pdb.set_trace()
            distDiagonal = D[i - 1][j - 1] + (s1[i - 1] != s2[j - 1])
            D[i][j] = min(distLeft, distAbove, distDiagonal)

    # return the last value (i.e., right most bottom value)
    return D[-1][-1]
