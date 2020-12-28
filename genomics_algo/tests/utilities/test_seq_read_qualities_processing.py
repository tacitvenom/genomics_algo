from genomics_algo.utilities.seq_read_qualities_processing import (
    map_phred33_to_error_probability,
    map_errorprobability_to_phred33,
    get_freq_for_qualities,
)


DELTA = 10e-4


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
