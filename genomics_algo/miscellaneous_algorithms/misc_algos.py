from itertools import product
import numpy as np

from typing import List, Set, Tuple

from genomics_algo.utilities.misc_utilities import (
    get_frequency_map,
    validate_bases_in_genome,
)
from genomics_algo.utilities.string_cmp import find_hamming_distance


def find_most_freq_k_substring(
    text: str, substring_length: int
) -> Tuple[List[str], int]:
    """
    Find the most frequent substring of length in a given text
    >>> find_most_freq_k_substring("GTACGTACC", 1)
    (['C'], 3)
    >>> find_most_freq_k_substring("GTACGTACC", 2)
    (['GT', 'TA', 'AC'], 2)
    >>> find_most_freq_k_substring("GTACGTACC", 4)
    (['GTAC'], 2)
    >>> find_most_freq_k_substring("GTACGTACC", 6)
    (['GTACGT', 'TACGTA', 'ACGTAC', 'CGTACC'], 1)
    """
    freq_map = get_frequency_map(text=text, substring_length=substring_length)
    frequency = max(freq_map.values())
    frequent_substrings = [key for key, value in freq_map.items() if value == frequency]
    return frequent_substrings, frequency


def find_pattern_clumps(
    text: str, substring_length: int, window_length: int, minimum_frequency: int
) -> Set[str]:
    """Find patterns forming clumps in a `text`, i.e., returns all the substrings of
    length `substring_length` in `text` which occurred at least `minimum_frequency` times
    in a window of fixed length `window_length` along the `text`, essentially looking for
    a region where a k-mer appears several times in short succession

    Returns:
        Set[str]: set of strings
    """
    patterns = set()
    for index in range(len(text) - window_length + 1):
        window = text[index : index + window_length]
        freq_map = get_frequency_map(text=window, substring_length=substring_length)
        for key, value in freq_map.items():
            if value >= minimum_frequency:
                patterns.add(key)
    return patterns


def find_minimum_gc_skew_location(genome: str) -> int:
    """TODO: [summary]

    Returns:
        [type]: [description]
    """
    assert set(genome) - {"A", "C", "G", "T"} == set()
    gc_skew = np.zeros(len(genome) + 1)
    for index in range(len(genome)):
        gc_skew[index + 1] = gc_skew[index]
        gc_skew[index + 1] += genome[index] == "G"
        gc_skew[index + 1] -= genome[index] == "C"
    return np.where(gc_skew == gc_skew.min())[0] - 1


def find_frequent_kmers_with_mismatches(genome: str, k: int, d: int) -> Set[str]:
    """Determine most frequent k-mers with at most `d` mismatches.
    A most frequent k-mer with up to `d` mismatches in `genome` is simply a string pattern maximising
    the total number of occurrences of said pattern in `genome` with at most `d` mismatches.
    Note that the pattern does not need to actually appear as a substring of `genome`.
    >>> find_frequent_kmers_with_mismatches('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)-{'ATGC', 'GATG', 'ATGT'}
    set()

    Parameters
    ----------
    genome : str
        String representation of genome.
    k: int
        Length of kmers to find.
    d: int
        Number of allowed mismatches in kmers.

    Returns
    -------
    Set[str]
        Set of most frequent kmers with up to d mismatches
    """

    n = len(genome)
    chars = {"A", "C", "G", "T"}
    # input validation:
    validate_bases_in_genome(genome)
    if n < k or k < d or d < 0:
        raise ValueError(
            f"The input values for genome, k and d don't make sense. It must hold: len(genome)>=k, k>=d, d>=0. Received: len(genome)={n}, k={k}, d={d}."
        )
    if k > 12 or d > 3:
        raise Warning(
            f"The large input values k={k} and/or d={d} might cause long run times."
        )

    frequency_map = {}
    # FIXME here ALL possible patterns of length k are created -> should be optimised
    possible_patterns = ["".join(p) for p in product(chars, repeat=k)]
    for i in range(n - k + 1):
        pattern = genome[i : i + k]
        for kmer in possible_patterns:
            if find_hamming_distance(pattern, kmer) <= d:
                if kmer in frequency_map.keys():
                    frequency_map[kmer] += 1
                else:
                    frequency_map[kmer] = 1

    most_frequent = max(frequency_map.values())
    return {
        kmer for kmer, frequency in frequency_map.items() if frequency == most_frequent
    }
