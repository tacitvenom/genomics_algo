import numpy as np

from typing import Dict, List, Tuple

from genomics_algo.helper import (
    get_frequency_map,
    reverse_complement,
    longest_common_suffix,
)


def get_occurences_with_naive_match(pattern: str, text: str) -> List[int]:
    """Get indices of all occurences of the string `pattern` in the
    string `text` using naive matching
    """
    occurences = []
    len_pattern = len(pattern)
    len_text = len(text)
    if len_pattern <= len_text:
        for index in range(len_text - len_pattern + 1):
            match = all(
                pattern[offset] == text[index + offset] for offset in range(len_pattern)
            )

            if match:
                occurences.append(index)

    return occurences


def get_occurences_with_exact_match_with_reverse_complement(
    pattern: str, text: str, exact_matching_algo: callable
) -> List[int]:
    """Get indices of all occurences of the DNA strand string `pattern` in the
    string `text` using naive matching considering reverse complement of the
    `pattern` string also
    """
    occurences = []
    occurences += exact_matching_algo(pattern, text)
    pattern_reverse_complement = reverse_complement(pattern)
    if pattern_reverse_complement != pattern:
        occurences += exact_matching_algo(pattern_reverse_complement, text)
    return occurences


def _get_alignments_skipped_bad_char_rule(
    mismatched_char: str, pattern_prefix: str
) -> int:
    """Get the number of alignments that can be skipped according to bad
    character rule in Boyer Moore's exact matching algorithm
    >>> _get_alignments_skipped_bad_char_rule("C", "ATCTTTATCATA")
    3
    >>> _get_alignments_skipped_bad_char_rule("G", "ATCTTTATCATA")
    12
    >>> _get_alignments_skipped_bad_char_rule("T", "GTAGCGGC")
    6
    >>> _get_alignments_skipped_bad_char_rule("C", "GTAGC")
    0
    >>> _get_alignments_skipped_bad_char_rule("C", "GT")
    2
    """
    reverse_pattern_prefix = pattern_prefix[::-1]  # checking from right to left
    len_pattern_prefix = len(pattern_prefix)
    for index in range(len_pattern_prefix):
        if reverse_pattern_prefix[index] == mismatched_char:
            return index
    return len(pattern_prefix)


def _get_alignments_skipped_good_suffix_rule(matched_suffix: str, pattern: str) -> int:
    """Get the number of alignments that can be skipped according to good
    suffix rule in Boyer Moore's exact matching algorithm
    >>> _get_alignments_skipped_good_suffix_rule("", "GTAGCGGCG")
    0
    >>> _get_alignments_skipped_good_suffix_rule("GCG", "GTAGCGGCG")
    2
    >>> _get_alignments_skipped_good_suffix_rule("GCGGCG", "GTAGCGGCG")
    7
    >>> _get_alignments_skipped_good_suffix_rule("TAC", "CTTACTTAC")
    3
    >>> _get_alignments_skipped_good_suffix_rule("TACTTAC", "CTTACTTAC")
    3
    """
    len_pattern = len(pattern)
    for i in range(len_pattern - 1, 0, -1):
        len_smaller_substr = min(i, len(matched_suffix))
        if (
            longest_common_suffix(pattern[:i], matched_suffix)
            == matched_suffix[-len_smaller_substr::]
        ):
            return len_pattern - i - 1
    return len_pattern - 1


def _get_alignments_skipped_gs_lookup(pattern: str) -> Dict[str, int]:
    """Get the number of alignments that can be skipped according to good
    suffix rule in Boyer Moore's exact matching algorithm for each possible
    suffix of a pattern in a dictionary
    """
    gs_lookup = {}
    for index in range(len(pattern)):
        suffix = pattern[index + 1 :]
        gs_lookup[suffix] = _get_alignments_skipped_good_suffix_rule(
            matched_suffix=suffix, pattern=pattern
        )
    return gs_lookup


def _get_alignments_skipped_bc_lookup(pattern: str) -> Dict[str, Dict[str, int]]:
    """Get the number of alignments that can be skipped according to bad
    character rule in Boyer Moore's exact matching algorithm for each possible
    character of the string with every possible prefix of a pattern in a dictionary
    of dictionaries
    """
    bc_lookup = {}
    unique_characters_in_pattern = set(pattern)
    for character in unique_characters_in_pattern:
        char_lookup = {}
        for index in range(len(pattern)):
            prefix = pattern[:index]
            char_lookup[prefix] = _get_alignments_skipped_bad_char_rule(
                mismatched_char=character, pattern_prefix=prefix
            )
        bc_lookup[character] = char_lookup
    return bc_lookup


def get_occurences_with_boyer_moore_exact_matching(
    pattern: str, text: str
) -> List[int]:
    """Get indices of all occurences of the string `pattern` in the
    string `text` using boyer-moore's exact matching
    """
    occurences = []
    len_pattern = len(pattern)
    len_text = len(text)
    if len_pattern <= len_text:
        index = 0
        while index < (len_text - len_pattern + 1):
            match = True
            for offset in range(len_pattern - 1, -1, -1):
                if pattern[offset] != text[index + offset]:
                    match = False
                    alignments_to_skip_bc = _get_alignments_skipped_bad_char_rule(
                        mismatched_char=text[index + offset],
                        pattern_prefix=pattern[:offset],
                    )
                    alignments_to_skip_gs = _get_alignments_skipped_good_suffix_rule(
                        matched_suffix=pattern[offset + 1 :], pattern=pattern
                    )
                    alignments_to_skip = max(
                        alignments_to_skip_bc, alignments_to_skip_gs
                    )
                    index += alignments_to_skip
                    break
            if match:
                occurences.append(index)
            index += 1
    return occurences


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


# def find_pattern_clumps(
#     text: str, substring_length: int, window_length: int, minimum_frequency: int
# ):
#     patterns = set()
#     for index in range(len(text) - window_length + 1):
#         window = text[index : index + window_length]
#         freq_map = get_frequency_map(text=window, substring_length=substring_length)
#         for key, value in freq_map.items():
#             if value >= minimum_frequency:
#                 patterns.add(key)
#     return patterns


def find_pattern_clumps(
    text: str, substring_length: int, window_length: int, minimum_frequency: int
):
    window = text[:window_length]
    initial_substring = window[:substring_length]
    freq_map = get_frequency_map(text=window, substring_length=substring_length)
    patterns = {
        key for key, value in freq_map.items() if value >= minimum_frequency
    }

    print("0", freq_map)
    for index in range(1, len(text) - window_length + 1):
        window = text[index : index + window_length]
        previous_freq_map = freq_map.copy()
        freq_map = {}
        last_substring = window[::-1][:substring_length][
            ::-1
        ]  # last `substring_length` characters
        for key, value in previous_freq_map.items():
            if key == initial_substring:
                freq_map[key] = previous_freq_map[key] - 1
            elif key == last_substring:
                freq_map[key] = previous_freq_map[key] + 1
            else:
                freq_map[key] = previous_freq_map[key]
        if last_substring not in freq_map:
            freq_map[last_substring] = 1
        for key, value in freq_map.items():
            if value >= minimum_frequency:
                patterns.add(key)
        initial_substring = window[:substring_length]
        print(index, freq_map)
    return patterns


def find_minimum_gc_skew_location(genome: str) -> int:
    assert set(genome) - {"A", "C", "G", "T"} == set()
    gc_skew = np.zeros(len(genome) + 1)
    for index in range(len(genome)):
        gc_skew[index + 1] = gc_skew[index]
        gc_skew[index + 1] += genome[index] == "G"
        gc_skew[index + 1] -= genome[index] == "C"
    return np.where(gc_skew == gc_skew.min())[0] - 1
