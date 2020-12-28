import numpy as np

from typing import List, Tuple

from genomics_algo.utilities.misc_utilities import get_frequency_map


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
):
    patterns = set()
    for index in range(len(text) - window_length + 1):
        window = text[index : index + window_length]
        freq_map = get_frequency_map(text=window, substring_length=substring_length)
        for key, value in freq_map.items():
            if value >= minimum_frequency:
                patterns.add(key)
    return patterns


def find_minimum_gc_skew_location(genome: str) -> int:
    assert set(genome) - {"A", "C", "G", "T"} == set()
    gc_skew = np.zeros(len(genome) + 1)
    for index in range(len(genome)):
        gc_skew[index + 1] = gc_skew[index]
        gc_skew[index + 1] += genome[index] == "G"
        gc_skew[index + 1] -= genome[index] == "C"
    return np.where(gc_skew == gc_skew.min())[0] - 1
