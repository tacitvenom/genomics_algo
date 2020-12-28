import math

from collections import Counter
from typing import List, Tuple


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
