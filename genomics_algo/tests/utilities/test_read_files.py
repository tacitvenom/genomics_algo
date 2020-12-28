from genomics_algo.utilities.read_files import (
    read_genome,
    read_fastq,
)


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
