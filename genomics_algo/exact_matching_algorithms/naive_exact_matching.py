from typing import List

from genomics_algo.utilities.misc_utilities import reverse_complement


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
