import math

COMPLEMENTARY_BASE = {"A": "T", "C": "G", "G": "C", "T": "A"}


def longest_common_prefix(s1, s2):
    """
    Finds the longest common prefix (substring) given two strings

    s1: string
        First string to compare
    s2: string
        Second string to compare

    Returns:
    longest_common_prefix: String
        Longest common prefix between s1 and s2
    """
    i = 0
    while i < min(len(s1), len(s2)):
        if s1[i] != s2[i]:
            break
        i += 1
    return s1[:i]


def reverse_complement(s):
    """
    Find the reverse complement of a DNA strand
    s: string
        A DNA sequence of a strand - the string must have the characters A, C, G and T

    Returns:
        DNA sequence of the opposite strand in the reverse order
    """
    # TODO: make robust againt garbage values
    return "".join([COMPLEMENTARY_BASE[base] for base in s[::-1]])


def read_genome(filename):
    """
    Reads a genome from a .fa file

    filename: string
        relative or absolute path of the .fa file to be read from

    Returns:
    genome: string
        Genome string
    """
    with open(filename) as f:
        genome = "".join(
            [line.rstrip() for line in f.readlines() if not line.startswith(">")]
        )
    return genome


def read_fastq(filename):
    """
    Reads sequences and qualities from a .fastq file

    filename: string
        relative or absolute path of the .fa file to be read from

    Returns:
    reads: list
        List of sequence reads
    qualities: list
        List of qualities corresponding to each sequence read

    """
    reads = []
    qualities = []
    with open(filename, "r") as f:
        while True:
            f.readline()
            read = f.readline().rstrip()
            f.readline()
            seq_qualities = f.readline().rstrip()
            if len(read) == 0:
                break
            reads.append(read)
            qualities.append(seq_qualities)

    return reads, qualities


def map_phred33_to_error_probability(phred33):
    """Maps a ASCII phred33 quality character to a quality score in fraction
    # TODO
    Args:
        quality ([type]): [description]
    """
    return 10 ** (-phred33 / 10)


def map_errorprobability_to_phred33(probability):
    """
    # TODO
    """
    return -10 * math.log(probability, 10)


def map_phred33_ascii_to_qualityscore(phred33_char):
    """Maps a ASCII phred33 quality character to a quality score
    # TODO
    Args:
        quality ([type]): [description]
    """
    return ord(phred33_char) - 33
