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
            distLeft = D[i][j - 1] + 1  # deletion in pattern
            distAbove = D[i - 1][j] + 1  # insertion in pattern
            distDiagonal = D[i - 1][j - 1] + (s1[i - 1] != s2[j - 1])  # substitution
            D[i][j] = min(distLeft, distAbove, distDiagonal)

    # return the last value (i.e., right most bottom value)
    return D[-1][-1]
