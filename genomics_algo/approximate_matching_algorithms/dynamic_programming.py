import doctest

from typing import List


def get_occurences_with_dynamic_programming(
    pattern: str, text: str, max_mismatches: int
) -> List[int]:
    """Get indices of all occurences of the string `pattern` in the string `text` using
    approximate matching (Levenshtein distance is used to count the number of mismatches i.e.,
    minimum number of edits including substitution, insertion and deletion needed in a string to
    turn it into another)
    >>> get_occurences_with_dynamic_programming("GCGTATGC", "TATTGGCTATACGGTT", 1) # doctest: +SKIP
    []
    >>> get_occurences_with_dynamic_programming("GCGTATGC", "TATTGGCTATACGGTT", 2) # doctest: +SKIP
    [5]
    >>> get_occurences_with_dynamic_programming("ACT", "GACTACGGAGACT", 0) # doctest: +SKIP
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
            distLeft = D[i][j - 1] + 1  # deletion in pattern
            distAbove = D[i - 1][j] + 1  # insertion in pattern
            distDiagonal = D[i - 1][j - 1] + (
                pattern[i - 1] != text[j - 1]
            )  # substitution
            D[i][j] = min(distLeft, distAbove, distDiagonal)

    # minimum in the bottom-most row should be at most the value of `max_mismatches`
    if min(D[-1]) > max_mismatches:
        return []
    else:
        occurence_end_indices = []
        for end_index, mismatch_count in enumerate(D[-1]):
            if mismatch_count <= max_mismatches:
                occurence_end_indices.append(end_index)

        occurences = _backtrace_approximate_match(D, occurence_end_indices)

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
    raise NotImplementedError("Helper not implemented yet!")
