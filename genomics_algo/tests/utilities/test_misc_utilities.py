import pytest

from genomics_algo.utilities.misc_utilities import (
    reverse_complement,
    generate_artificial_reads,
    get_frequency_map,
)


def test_reverse_complement():
    with pytest.raises(KeyError):
        reverse_complement("AGRTTTAG")


def test_generate_artificial_reads():
    genome = "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC"
    read_length = 5
    number_of_reads = 3

    reads = generate_artificial_reads(
        genome=genome, number_of_reads=number_of_reads, read_length=read_length
    )

    for read in reads:
        assert read in genome


def test_get_frequency_map():
    with pytest.raises(AssertionError):
        get_frequency_map("", 1)
    with pytest.raises(AssertionError):
        get_frequency_map("GTACGTACC", 0)
    with pytest.raises(AssertionError):
        get_frequency_map("GTACGTACC", -2)
