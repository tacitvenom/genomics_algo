from typing import List, Tuple


def read_genome(filename: str) -> str:
    """
    Reads a genome from a .fa file

    filename: relative or absolute path of the .fa file to be read from

    Returns:
        Genome string
    """
    with open(filename) as f:
        genome = "".join(
            [line.rstrip() for line in f.readlines() if not line.startswith(">")]
        )
    return genome


def read_fastq(filename: str) -> Tuple[List[str], List[str]]:
    """
    Reads sequences and qualities from a .fastq file

    filename: relative or absolute path of the .fa file to be read from

    Returns:
        List of sequence reads
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
