from genomics_algo.algorithms import get_occurences_with_naive_match


def test_get_occurences_with_naive_match():
    text = "GACTACGGAGACT"
    pattern = "ACT"

    result = get_occurences_with_naive_match(pattern, text)
    assert result == [1, 10]
