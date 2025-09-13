from enum import Enum
from typing import Dict, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SeqOrientation(int, Enum):
    """
    Orientation of a sequence (forward or reverse compliment) relative to some reference
    sequence.
    """

    FORWARD = 1
    REVERSE = -1

    def __mul__(self, other: "SeqOrientation | int") -> "SeqOrientation":
        if not isinstance(other, SeqOrientation):
            other = SeqOrientation(other)
        return SeqOrientation(self.value * other.value)

    def __rmul__(self, other: "SeqOrientation | int") -> "SeqOrientation":
        if not isinstance(other, SeqOrientation):
            other = SeqOrientation(other)
        return SeqOrientation(self.value * other.value)


def is_subseq(needle: str | Seq, haystack: str | Seq, check_rc: bool = True) -> bool:
    """Check whether `needle` is a substring of `haystack` without errors.

    Args:
        needle: The sequence to search for, as either a string or Seq object
        haystack: The sequence to search within, as either a string or Seq object
        check_rc: If True, check both the forward sequence and its reverse complement

    Returns:
        bool: True if needle (or its reverse complement when check_rc=True) is found
              within haystack
    """
    # string operations are faster than Seq operations
    needle_str = str(needle)
    haystack_str = str(haystack)

    if not check_rc:
        return needle_str in haystack_str

    return (
        needle_str in haystack_str
        or str(Seq(needle_str).reverse_complement()) in haystack_str
    )


def contig_ids_by_seed(
    records: List[SeqRecord], seed_seqs: List[Seq]
) -> Dict[int, SeqOrientation]:
    """
    Given a list of contigs and a list of seed sequences, returns the indices of contigs
    that contain at least one seed, along with the orientation of the contig with respect
    to the seed.

    Args:
        records: List of SeqRecord objects representing contigs
        seed_seqs: List of seed sequences to search for
    Returns:
        Dict whose keys correspond to the indices of records that contain seed sequences,
        and whose values correspond to the contig orientation with respect to the seed
    """
    filtered_records = {}
    for i, rec in enumerate(records):
        contig_sequence = str(rec.seq)
        for seed in seed_seqs:
            if is_subseq(needle=str(seed), haystack=contig_sequence, check_rc=False):
                filtered_records[i] = SeqOrientation.FORWARD
                break
            elif is_subseq(
                needle=str(seed.reverse_complement()),
                haystack=contig_sequence,
                check_rc=False,
            ):
                filtered_records[i] = SeqOrientation.REVERSE
                break

    return filtered_records
