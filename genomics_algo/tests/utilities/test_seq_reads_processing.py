import pytest

import numpy as np


from genomics_algo.utilities.seq_reads_processing import (
    same_length_reads,
    find_GC_by_position,
)


def test_same_length_reads():
    with pytest.raises(AssertionError):
        same_length_reads([])


def test_find_GC_by_position():
    reads = []
    with pytest.raises(AssertionError):
        find_GC_by_position(reads)

    with pytest.raises(AssertionError):
        find_GC_by_position(["ACGTNO"])

    reads = ["AACGTTA", "CGCGTTT", "GTTAC"]
    with pytest.raises(AssertionError):
        find_GC_by_position(reads)

    reads = ["AACGTTA"]
    assert all(find_GC_by_position(reads) == np.array([0, 0, 1, 1, 0, 0, 0]))

    reads = ["AACGTTA", "CGCGTTT"]
    assert all(find_GC_by_position(reads) == np.array([0.5, 0.5, 1, 1, 0, 0, 0]))
