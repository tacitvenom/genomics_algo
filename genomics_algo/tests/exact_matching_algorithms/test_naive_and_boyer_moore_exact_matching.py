import pytest

from genomics_algo.utilities.read_files import read_genome

from genomics_algo.exact_matching_algorithms.boyer_moore_exact_matching import (
    get_occurences_with_boyer_moore_exact_matching,
)
from genomics_algo.exact_matching_algorithms.naive_exact_matching import (
    get_occurences_with_naive_match,
    get_occurences_with_exact_match_with_reverse_complement,
)


@pytest.mark.parametrize(
    "exact_matching_algo",
    [get_occurences_with_naive_match, get_occurences_with_boyer_moore_exact_matching],
)
def test_get_occurences_with_exact_match(exact_matching_algo):
    text = "GACTACGGAGACT"
    pattern = "ACT"
    result = exact_matching_algo(pattern, text)
    assert result == [1, 10]

    text = "there would have been a time for such a word"
    pattern = "word"
    result = exact_matching_algo(pattern, text)
    assert result == [40]

    text = "needle need noodle needle"
    pattern = "needle"
    result = exact_matching_algo(pattern, text)
    assert result == [0, 19]

    text = "GCTTCTGCTACCTTATGCGCGCGCCTTTTGCCGCGGACCTTTTGCA"
    pattern = "CCTTTTGC"
    result = exact_matching_algo(pattern, text)
    assert result == [23, 37]

    text = "CGATATATCCATAG"
    pattern = "ATA"
    result = exact_matching_algo(pattern, text)
    assert result == [2, 4, 10]


@pytest.mark.skip(reason="Overhead of 2-3 seconds")
@pytest.mark.parametrize(
    "exact_matching_algo",
    [get_occurences_with_naive_match, get_occurences_with_boyer_moore_exact_matching],
)
def test_get_occurences_in_entire_genome_with_boyer_moores_exact_match(
    exact_matching_algo,
):
    text = read_genome("genomics_algo/tests/test_data/vibrio_cholerae.txt")
    pattern = "ATGATCAAG"
    result = exact_matching_algo(pattern, text)
    assert result == [
        116556,
        149355,
        151913,
        152013,
        152394,
        186189,
        194276,
        200076,
        224527,
        307692,
        479770,
        610980,
        653338,
        679985,
        768828,
        878903,
        985368,
    ]


@pytest.mark.parametrize(
    "exact_matching_algo",
    [get_occurences_with_naive_match, get_occurences_with_boyer_moore_exact_matching],
)
def test_get_occurences_with_exact_match_with_reverse_complement(exact_matching_algo):
    text = "GACTACGGAGACT"
    pattern = "ACT"
    result = get_occurences_with_exact_match_with_reverse_complement(
        pattern, text, exact_matching_algo
    )
    assert result == [1, 10]

    text = "AAAAAAAAAACCCAAAAAAAAAAGGGAAAAAAAAAA"
    pattern = "CCC"
    result = get_occurences_with_exact_match_with_reverse_complement(
        pattern, text, exact_matching_algo
    )
    assert result == [10, 23]

    text = "AAAAAAAAAACGCGAAAAAAAAAACGCGAAAAAAAAAA"
    pattern = "CGCG"
    result = get_occurences_with_exact_match_with_reverse_complement(
        pattern, text, exact_matching_algo
    )
    assert result == [10, 24]

    text = read_genome("genomics_algo/tests/test_data/phix.fa")
    pattern = "ATTA"
    result = get_occurences_with_exact_match_with_reverse_complement(
        pattern, text, exact_matching_algo
    )
    assert min(result) == 62
    assert len(result) == 60
