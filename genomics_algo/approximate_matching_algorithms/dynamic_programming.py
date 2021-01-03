import doctest

from typing import List


def get_occurences_with_dynamic_programming(
    pattern: str, text: str, max_mismatches: int
) -> List[int]:
    """Get indices of all occurences of the string `pattern` in the string `text` using
    approximate matching (Levenshtein distance is used to count the number of mismatches i.e.,
    minimum number of edits including substitution, insertion and deletion needed in a string to
    turn it into another)
    >>> get_occurences_with_dynamic_programming("GCGTATGC", "TATTGGCTATACGGTT", 1)
    []
    >>> get_occurences_with_dynamic_programming("GCGTATGC", "TATTGGCTATACGGTT", 2)
    [5]
    >>> get_occurences_with_dynamic_programming("ACT", "GACTACGGAGACT", 0)
    [1, 10]
    """
    assert len(pattern) <= len(text)

    # initializing a matrix for with `len(pattern) + 1` rows and `len(text) + 1` columns
    D = [[0 for x in range(len(text) + 1)] for y in range(len(pattern) + 1)]

    # fill first column
    for i in range(len(pattern) + 1):
        D[i][0] = i
    # fill first row
    for j in range(len(text) + 1):
        D[0][j] = 0

    # fill rest of the matrix
    for i in range(1, len(pattern) + 1):
        for j in range(1, len(text) + 1):
            distance_left = D[i][j - 1] + 1  # deletion in pattern
            distance_above = D[i - 1][j] + 1  # insertion in pattern
            distance_diagonal = D[i - 1][j - 1] + (
                pattern[i - 1] != text[j - 1]
            )  # substitution
            D[i][j] = min(distance_left, distance_above, distance_diagonal)

    # minimum in the bottom-most row should be at most the value of `max_mismatches`
    if min(D[-1]) > max_mismatches:
        return []
    else:
        occurence_end_indices = []
        for end_index, mismatch_count in enumerate(D[-1]):
            if mismatch_count <= max_mismatches:
                occurence_end_indices.append(end_index)

        occurences = _backtrace_approximate_match(
            pattern=pattern, text=text, D=D, occurence_end_indices=occurence_end_indices
        )

    return occurences


def _backtrace_approximate_match(
    pattern: str, text: str, D: List[List[int]], occurence_end_indices: List[int]
) -> List[int]:
    """Helper function that backtraces and find the beginning index of an approximate occurence
    provided the dynamic programmically filled matrix `D` and the list of end indices of an approximate
    match

    Args:
        pattern (str): `pattern` to be searched in the `text`
        text (str): `text` in which `pattern` is searched
        D (List[List[int]]): Matrix with `len(pattern) + 1` rows and `len(text) + 1` columns with mismatch
            information using dynamic programming
        occurence_end_indices (List[int]): List of indices of the end of the occurence of an approximate
            match of `pattern` in the `text`

    Returns:
        List[int]: List of integers where the `pattern` approximately matches in the `text`
    """
    occurence_start_indices = []
    for occurence_end_index in occurence_end_indices:
        i = len(D) - 1
        j = occurence_end_index
        while i > 0:
            distance_value = D[i][j]
            distance_left = D[i][j - 1]
            distance_above = D[i - 1][j]
            distance_diagonal = D[i - 1][j - 1]
            if distance_value == distance_above + 1:
                i -= 1
            elif distance_value == distance_left + 1:
                j -= 1
            elif distance_value == distance_diagonal + (pattern[i - 1] != text[j - 1]):
                i -= 1
                j -= 1
        occurence_start_indices.append(j)
    return occurence_start_indices
