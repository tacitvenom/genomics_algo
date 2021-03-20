import numpy as np
import pytest

from genomics_algo.utilities.read_files import read_genome
from genomics_algo.miscellaneous_algorithms.misc_algos import (
    find_pattern_clumps,
    find_minimum_gc_skew_location,
    find_frequent_kmers_with_mismatches,
)


@pytest.mark.skip(reason="Takes 10-15 mins in current implementation")
def test_find_pattern_clumps_with_genome():
    text = read_genome("genomics_algo/tests/test_data/genomes/e_coli.txt")
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
    genome = read_genome("genomics_algo/tests/test_data/genomes/e_coli.txt")
    result = find_minimum_gc_skew_location(genome)
    np.testing.assert_array_equal([3923619, 3923620, 3923621, 3923622], result)


def test_find_frequent_kmers_with_mismatches_raises():
    with pytest.raises(Warning):
        find_frequent_kmers_with_mismatches("ACGTTGCAACGTTGCA", 13, 3)
    with pytest.raises(Warning):
        find_frequent_kmers_with_mismatches("ACGTTGCAACGTTGCA", 12, 4)
    with pytest.raises(ValueError) as e:
        find_frequent_kmers_with_mismatches("ACGT", 6, -1)
        assert "Received: len(genome)=4, k=6, d=-1." in str(e.value)


@pytest.mark.skip(reason="Takes around 2 minutes to execute")
def test_find_frequent_kmers_with_mismatches_benchmark():
    res = find_frequent_kmers_with_mismatches("ACGTTGCAACGTTGCA", 12, 3)
    assert len(res) == 32855


def test_find_frequent_kmers_with_mismatches():
    """Some debug datasets taken from:
    http://bioinformaticsalgorithms.com/data/debugdatasets/replication/FrequentWordsWithMismatchesProblem.pdf
    """

    """ Dataset 1
    This dataset checks that the implementation includes k-mers that do not actually appear in Text.
    Notice here that, although none of the output k-mers except for AA actually appear in Text, they
    are all valid because they appear in Text with up to 1 mismatch (i.e. 0 or 1 mismatch).
    """
    genome1 = "AAAAAAAAAA"
    k1 = 2
    d1 = 1
    expected1 = {"AA", "AC", "AG", "CA", "AT", "GA", "TA"}
    result1 = find_frequent_kmers_with_mismatches(genome1, k1, d1)
    assert result1 == expected1

    """ Dataset 2
    This dataset makes sure that the implementation is not accidentally swapping k and d.
    """
    genome2 = "AGTCAGTC"
    k2 = 4
    d2 = 2
    expected2 = {
        "TCTC",
        "CGGC",
        "AAGC",
        "TGTG",
        "GGCC",
        "AGGT",
        "ATCC",
        "ACTG",
        "ACAC",
        "AGAG",
        "ATTA",
        "TGAC",
        "AATT",
        "CGTT",
        "GTTC",
        "GGTA",
        "AGCA",
        "CATC",
    }
    result2 = find_frequent_kmers_with_mismatches(genome2, k2, d2)
    assert result2 == expected2

    """ Dataset 3
    This dataset makes sure you are not finding patterns in the Reverse Complement of genome
    """
    genome3 = "AATTAATTGGTAGGTAGGTA"
    k3 = 4
    d3 = 0
    expected3 = {"GGTA"}
    result3 = find_frequent_kmers_with_mismatches(genome3, k3, d3)
    assert result3 == expected3

    """ Dataset 4
    This dataset first checks that k-mers with exactly d mismatches are being found. Then, it
    checks that k-mers with less than d mismatches are being allowed (i.e. you are not only allowing
    k-mers with exactly d mismatches). Next, it checks that you are not returning too few k-mers.
    Last, it checks that you are not returning too many k-mers.
    """
    genome4 = "ATA"
    k4 = 3
    d4 = 1
    expected4 = {"GTA", "ACA", "AAA", "ATC", "ATA", "AGA", "ATT", "CTA", "TTA", "ATG"}
    result4 = find_frequent_kmers_with_mismatches(genome4, k4, d4)
    assert result4 == expected4

    """ Dataset 5
    This dataset checks that your code is not looking for k-mers in the Reverse Complement
    of genome.
    """
    genome5 = "AAT"
    k5 = 3
    d5 = 0
    expected5 = {"AAT"}
    result5 = find_frequent_kmers_with_mismatches(genome5, k5, d5)
    assert result5 == expected5
