from genomics_algo.helper import read_genome
from genomics_algo.algorithms import (
    get_occurences_with_naive_match,
    get_occurences_with_naive_match_with_reverse_complement,
)


def test_get_occurences_with_naive_match():
    text = "GACTACGGAGACT"
    pattern = "ACT"
    result = get_occurences_with_naive_match(pattern, text)
    assert result == [1, 10]


def test_get_occurences_with_naive_match_with_reverse_complement():
    text = "GACTACGGAGACT"
    pattern = "ACT"
    result = get_occurences_with_naive_match_with_reverse_complement(pattern, text)
    assert result == [1, 10]

    text = "AAAAAAAAAACCCAAAAAAAAAAGGGAAAAAAAAAA"
    pattern = "CCC"
    result = get_occurences_with_naive_match_with_reverse_complement(pattern, text)
    assert result == [10, 23]

    text = "AAAAAAAAAACGCGAAAAAAAAAACGCGAAAAAAAAAA"
    pattern = "CGCG"
    result = get_occurences_with_naive_match_with_reverse_complement(pattern, text)
    assert result == [10, 24]

    text = read_genome("genomics_algo/tests/test_data/phix.fa")
    pattern = "ATTA"
    result = get_occurences_with_naive_match_with_reverse_complement(pattern, text)
    assert min(result) == 62
    assert len(result) == 60
