import numpy as np
import pytest

from collections import Counter

from genomics_algo.helper import (
    longest_common_prefix,
    longest_common_suffix,
    reverse_complement,
    read_genome,
    read_fastq,
    map_phred33_to_error_probability,
    map_errorprobability_to_phred33,
    map_phred33_ascii_to_qualityscore,
    get_freq_for_qualities,
    same_length_reads,
    find_GC_by_position,
    get_base_freq,
    generate_artificial_reads,
)

DELTA = 10e-4


def test_longest_common_prefix():
    s1 = "ACTA"
    s2 = "GCCT"
    assert longest_common_prefix(s1, s2) == ""
    s1 = "ACTA"
    s2 = "ACT"
    assert longest_common_prefix(s1, s2) == "ACT"
    s1 = "ACTA"
    s2 = "ACT"
    assert longest_common_prefix(s1, s2) == "ACT"
    s1 = "GATA"
    s2 = "GAAT"
    assert longest_common_prefix(s1, s2) == "GA"
    s1 = "ATGA"
    s2 = ""
    assert longest_common_prefix(s1, s2) == ""
    s1 = ""
    s2 = "GCCT"
    assert longest_common_prefix(s1, s2) == ""
    s1 = "GCCT"
    s2 = s1
    assert longest_common_prefix(s1, s2) == s1


def test_longest_common_suffix():
    s1 = "ACTA"
    s2 = "GCCT"
    assert longest_common_suffix(s1, s2) == ""
    s1 = "ACTA"
    s2 = "CTA"
    assert longest_common_suffix(s1, s2) == "CTA"
    s1 = "CTA"
    s2 = "ACTA"
    assert longest_common_suffix(s1, s2) == "CTA"
    s1 = "GATAT"
    s2 = "GAATAT"
    assert longest_common_suffix(s1, s2) == "ATAT"
    s1 = "ATGA"
    s2 = ""
    assert longest_common_suffix(s1, s2) == ""
    s1 = ""
    s2 = "GCCT"
    assert longest_common_suffix(s1, s2) == ""
    s1 = "GCCT"
    s2 = s1
    assert longest_common_prefix(s1, s2) == s1


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("") == ""
    with pytest.raises(KeyError):
        reverse_complement("AGRTTTAG")


def test_read_genome():
    genome = read_genome("genomics_algo/tests/test_data/lambda_virus.fa")
    assert len(genome) == 48502
    assert genome[:50] == "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAA"


def test_read_fastq():
    reads, qualities = read_fastq(
        "genomics_algo/tests/test_data/SRR835775_1.first1000.fastq"
    )
    assert len(reads) == len(qualities) == 1000
    assert (
        reads[0]
        == "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTCACCCTAACCCTAACCCTAACCGTATCCGTCACCCTAACCCTAAC"
    )
    assert (
        qualities[0]
        == "???B1ADDD8??BB+C?B+:AA883CEE8?C3@DDD3)?D2;DC?8?=BAD=@C@(.6.6=A?=?@##################################"
    )


def test_map_phred33_to_error_probability():
    phred33_to_error_probability = {
        10: 0.1,
        20: 0.01,
        50: 0.00001,
    }
    for score, probability in phred33_to_error_probability.items():
        assert abs(map_phred33_to_error_probability(score) - probability) < DELTA


def test_map_errorprobability_to_phred33():
    error_probability_to_phred33 = {
        0.1: 10,
        0.01: 20,
        0.00001: 50,
    }
    for probability, score in error_probability_to_phred33.items():
        assert abs(map_errorprobability_to_phred33(probability) - score) < DELTA


def test_map_phred33_ascii_to_qualityscore():
    assert map_phred33_ascii_to_qualityscore("#") == 2
    assert map_phred33_ascii_to_qualityscore("J") == 41


def test_get_freq_for_qualities():
    qualities = [
        "???B1ADDD8??BB+C?B+:AA883CEE8?C3@DDD3)",
        "JJJFJJFGIIIIH=CBFCF=CCEG)=>EHB2@@DEC18",
    ]
    scores, frequencies = get_freq_for_qualities(qualities=qualities)
    assert scores == [
        8,
        10,
        16,
        17,
        18,
        23,
        25,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
    ]
    assert frequencies == [
        2,
        2,
        2,
        1,
        3,
        5,
        1,
        3,
        1,
        7,
        3,
        3,
        6,
        8,
        7,
        5,
        4,
        2,
        2,
        4,
        5,
    ]


def test_same_length_reads():
    reads = []
    with pytest.raises(AssertionError):
        same_length_reads(reads)

    reads = ["AACGTTA"]
    assert same_length_reads(reads)

    reads = ["AACGTTA", "CGCGTTT"]
    assert same_length_reads(reads)

    reads = ["AACGTTA", "CGCGTTT", "GTTAC"]
    assert not same_length_reads(reads)


def test_find_GC_by_position():
    reads = []
    with pytest.raises(AssertionError):
        find_GC_by_position(reads)

    reads = ["AACGTTA", "CGCGTTT", "GTTAC"]
    with pytest.raises(AssertionError):
        find_GC_by_position(reads)

    reads = ["AACGTTA"]
    assert all(find_GC_by_position(reads) == np.array([0, 0, 1, 1, 0, 0, 0]))

    reads = ["AACGTTA", "CGCGTTT"]
    assert all(find_GC_by_position(reads) == np.array([0.5, 0.5, 1, 1, 0, 0, 0]))


def test_get_base_freq():
    reads = ["AACGTTAN"]
    assert get_base_freq(reads) == Counter({"G": 1, "C": 1, "A": 3, "T": 2, "N": 1})

    reads = ["AACGTTA", "CGCGTTT"]
    assert get_base_freq(reads) == Counter({"G": 3, "C": 3, "A": 3, "T": 5})


def test_generate_artificial_reads():
    genome = "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC"
    read_length = 5
    number_of_reads = 3

    reads = generate_artificial_reads(
        genome=genome, number_of_reads=number_of_reads, read_length=read_length
    )

    for read in reads:
        assert read in genome
