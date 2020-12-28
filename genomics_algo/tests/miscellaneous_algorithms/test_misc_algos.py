import numpy as np
import pytest

from genomics_algo.utilities.read_files import read_genome
from genomics_algo.miscellaneous_algorithms.misc_algos import (
    find_pattern_clumps,
    find_minimum_gc_skew_location,
)


@pytest.mark.skip(reason="Takes 10-15 mins in current implementation")
def test_find_pattern_clumps_with_genome():
    text = read_genome("genomics_algo/tests/test_data/e_coli.txt")
    patterns = find_pattern_clumps(
        text=text, substring_length=9, window_length=500, minimum_frequency=3
    )
    assert len(patterns) == 1904


def test_find_pattern_clumps_short():
    text = "GACAGAC"
    patterns = find_pattern_clumps(
        text=text, substring_length=3, window_length=7, minimum_frequency=2
    )
    assert patterns == {"GAC"}

    text = "GACCTACCGTATACGCCGACGACTTACTACATGCATGTAC"
    patterns = find_pattern_clumps(
        text=text, substring_length=3, window_length=16, minimum_frequency=3
    )
    assert patterns == {"TAC"}


@pytest.mark.skip(reason="Takes 20-25 seconds in current implementation")
def test_find_pattern_clumps_long():
    text = "GACCTACCGTATACGCCGACGACTTACTACATGCATGTAC" * 100_000
    patterns = find_pattern_clumps(
        text=text, substring_length=3, window_length=16, minimum_frequency=3
    )
    assert patterns == {"TAC"}


def test_find_minimum_gc_skew_location():
    genome = "CATGGGCATCGGCCATACGCCGAATA"
    assert find_minimum_gc_skew_location(genome) == 20
    genome = "CATGGGCATCGGCCATACGCCGAATACGA"
    result = find_minimum_gc_skew_location(genome)
    np.testing.assert_array_equal([20, 26], result)


@pytest.mark.skip(reason="Takes about 30s in current implementation")
def test_find_minimum_gc_skew_location_in_genome():
    genome = read_genome("genomics_algo/tests/test_data/e_coli.txt")
    result = find_minimum_gc_skew_location(genome)
    np.testing.assert_array_equal([3923619, 3923620, 3923621, 3923622], result)
