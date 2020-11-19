import pytest
from genomics_algo.helper import read_genome
from genomics_algo.algorithms import (
    get_occurences_with_naive_match,
    get_occurences_with_exact_match_with_reverse_complement,
    get_alignments_skipped_gs_lookup,
    get_occurences_with_boyer_moore_exact_matching,
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


def test_get_alignments_skipped_gs_lookup():
    pattern = ""
    assert get_alignments_skipped_gs_lookup(pattern) == {}
    pattern = "GTAGCGGCG"
    expected_lookup = {
        "": 0,
        "G": 1,
        "CG": 2,
        "GCG": 2,
        "GGCG": 7,
        "CGGCG": 7,
        "GCGGCG": 7,
        "AGCGGCG": 7,
        "TAGCGGCG": 7,
    }
    assert get_alignments_skipped_gs_lookup(pattern) == expected_lookup
