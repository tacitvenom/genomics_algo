import pytest

from genomics_algo.helper import (
    longest_common_prefix,
    reverse_complement,
    read_genome,
    read_fastq,
    map_phred33_to_error_probability,
    map_errorprobability_to_phred33,
    map_phred33_ascii_to_qualityscore,
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
