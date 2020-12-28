import numpy as np

from collections import Counter
from typing import List


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
