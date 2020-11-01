import pytest

from genomics_algo.helper import longest_common_prefix, reverse_complement, read_genome

def test_longest_common_prefix():
    s1 = 'ACTA'
    s2 = 'GCCT'
    assert longest_common_prefix(s1, s2) == ''
    s1 = 'ACTA'
    s2 = 'ACT'
    assert longest_common_prefix(s1, s2) == 'ACT'
    s1 = 'ACTA'
    s2 = 'ACT'
    assert longest_common_prefix(s1, s2) == 'ACT'
    s1 = 'GATA'
    s2 = 'GAAT'
    assert longest_common_prefix(s1, s2) == 'GA'
    s1 = 'ATGA'
    s2 = ''
    assert longest_common_prefix(s1, s2) == ''
    s1 = ''
    s2 = 'GCCT'
    assert longest_common_prefix(s1, s2) == ''


def test_reverse_complement():
    assert reverse_complement('ATGC') == 'GCAT'
    assert reverse_complement('') == ''
    with pytest.raises(KeyError):
        reverse_complement('AGRTTTAG')


def test_read_genome():
    genome = read_genome('genomics_algo/tests/test_data/lambda_virus.fa')
    assert len(genome) == 48502
