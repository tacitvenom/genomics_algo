from typing import Dict, List

from genomics_algo.utilities.string_cmp import longest_common_suffix


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
