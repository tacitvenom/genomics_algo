import numpy as np
import pytest

from genomics_algo.helper import read_genome
from genomics_algo.algorithms import (
    get_occurences_with_naive_match,
    get_occurences_with_exact_match_with_reverse_complement,
    _get_alignments_skipped_gs_lookup,
    _get_alignments_skipped_bc_lookup,
    get_occurences_with_boyer_moore_exact_matching,
    find_pattern_clumps,
    find_minimum_gc_skew_location,
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


def test__get_alignments_skipped_gs_lookup():
    pattern = ""
    assert _get_alignments_skipped_gs_lookup(pattern) == {}

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
    assert _get_alignments_skipped_gs_lookup(pattern) == expected_lookup


def test__get_alignments_skipped_bc_lookup():
    pattern = ""
    assert _get_alignments_skipped_bc_lookup(pattern=pattern) == {}

    pattern = "GTAGCGGCG"
    expected_lookup = {
        "A": {
            "": 0,
            "G": 1,
            "GT": 2,
            "GTA": 0,
            "GTAG": 1,
            "GTAGC": 2,
            "GTAGCG": 3,
            "GTAGCGG": 4,
            "GTAGCGGC": 5,
        },
        "C": {
            "": 0,
            "G": 1,
            "GT": 2,
            "GTA": 3,
            "GTAG": 4,
            "GTAGC": 0,
            "GTAGCG": 1,
            "GTAGCGG": 2,
            "GTAGCGGC": 0,
        },
        "G": {
            "": 0,
            "G": 0,
            "GT": 1,
            "GTA": 2,
            "GTAG": 0,
            "GTAGC": 1,
            "GTAGCG": 0,
            "GTAGCGG": 0,
            "GTAGCGGC": 1,
        },
        "T": {
            "": 0,
            "G": 1,
            "GT": 0,
            "GTA": 1,
            "GTAG": 2,
            "GTAGC": 3,
            "GTAGCG": 4,
            "GTAGCGG": 5,
            "GTAGCGGC": 6,
        },
    }
    assert _get_alignments_skipped_bc_lookup(pattern=pattern) == expected_lookup


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
