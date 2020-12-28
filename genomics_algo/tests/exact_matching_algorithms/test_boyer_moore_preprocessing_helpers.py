from genomics_algo.exact_matching_algorithms.boyer_moore_exact_matching import (
    _get_alignments_skipped_gs_lookup,
    _get_alignments_skipped_bc_lookup,
)


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
