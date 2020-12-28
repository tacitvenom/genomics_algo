import pytest

from genomics_algo.utilities.string_cmp import find_hamming_distance


def test_find_hamming_distance():
    with pytest.raises(AssertionError):
        find_hamming_distance("A", "ATG")
