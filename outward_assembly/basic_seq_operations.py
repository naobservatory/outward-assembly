from Bio.Seq import Seq


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
