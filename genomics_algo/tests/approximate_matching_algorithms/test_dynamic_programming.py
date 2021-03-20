from genomics_algo.approximate_matching_algorithms.dynamic_programming import (
    _backtrace_approximate_match,
)


def test__backtrace_approximate_match():
    text = "TATTGGCTATACGGTT"
    pattern = "GCGTATGC"
    D = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        [2, 2, 2, 2, 2, 1, 1, 0, 1, 2, 2, 2, 1, 1, 1, 1, 2],
        [3, 3, 3, 3, 3, 2, 1, 1, 1, 2, 3, 3, 2, 1, 1, 2, 2],
        [4, 3, 4, 3, 3, 3, 2, 2, 1, 2, 2, 3, 3, 2, 2, 1, 2],
        [5, 4, 3, 4, 4, 4, 3, 3, 2, 1, 2, 2, 3, 3, 3, 2, 2],
        [6, 5, 4, 3, 4, 5, 4, 4, 3, 2, 1, 2, 3, 4, 4, 3, 2],
        [7, 6, 5, 4, 4, 4, 5, 5, 4, 3, 2, 2, 3, 3, 4, 4, 3],
        [8, 7, 6, 5, 5, 5, 5, 5, 5, 4, 3, 3, 2, 3, 4, 5, 4],
    ]
    occurence_end_indices = [12]
    result_occurence_start_indices = _backtrace_approximate_match(
        pattern=pattern, text=text, D=D, occurence_end_indices=occurence_end_indices
    )
    expected_occurence_start_indices = [5]
    assert expected_occurence_start_indices == result_occurence_start_indices
