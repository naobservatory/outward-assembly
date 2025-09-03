import warnings
from typing import Dict, List, Sequence

import networkx as nx
from Bio.Seq import Seq


def _seqs_overlap_single(
    back_seq: str | Seq,
    front_seq: str | Seq,
    n_0_error: int,
    n_1_error: int,
    max_overlap: int | None = None,
) -> bool:
    """Return whether two sequences overlap by at least n_0_error bases with no errors
    or n_1_error bases with at most 1 error. Only considers overlaps where the end
    of the back sequence overlaps with the start of the forward sequence, like:

    back:    GATTACA
    forward:     ACATTG

    Note that no checks for reverse complements are done, nor does this function
    consider beginning of back_seq overlapping end of front_seq.

    Args:
        back_seq: The sequence that may overlap at its end
        front_seq: The sequence that may overlap at its start
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch
        max_overlap: Maximum overlap length to check (defaults to min length of sequences)

    Returns:
        True if sequences overlap with given criteria, False otherwise
    """
    back_seq = Seq(back_seq) if isinstance(back_seq, str) else back_seq
    front_seq = Seq(front_seq) if isinstance(front_seq, str) else front_seq

    if max_overlap is None:
        max_overlap = min(len(back_seq), len(front_seq))

    for overlap_len in range(n_0_error, max_overlap + 1):
        back_suffix = back_seq[-overlap_len:]
        front_prefix = front_seq[:overlap_len]
        disagreements = sum(b != f for b, f in zip(back_suffix, front_prefix))

        if disagreements == 0 or (overlap_len >= n_1_error and disagreements == 1):
            return True

    return False


def seqs_overlap(
    seq1: str | Seq, seq2: str | Seq, n_0_error: int, n_1_error: int
) -> bool:
    """Return whether two sequences overlap by at least n_0_error bases with no errors
    or n_1_error bases with at most 1 error. Detects reverse complement overlaps.

    Args:
        seq1: First sequence to check for overlap
        seq2: Second sequence to check for overlap
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch

    Returns:
        True if sequences overlap with given criteria, False otherwise

    Raises:
        ValueError: If n_1_error < n_0_error
    """
    if n_1_error < n_0_error:
        raise ValueError("n_1_error must be at least n_0_error")

    seq1 = Seq(seq1) if isinstance(seq1, str) else seq1
    seq2 = Seq(seq2) if isinstance(seq2, str) else seq2

    seq2_rc = seq2.reverse_complement()
    seq1_rc = seq1.reverse_complement()

    return (
        _seqs_overlap_single(seq1, seq2, n_0_error, n_1_error)
        or _seqs_overlap_single(seq1, seq2_rc, n_0_error, n_1_error)
        or _seqs_overlap_single(seq1_rc, seq2, n_0_error, n_1_error)
        or _seqs_overlap_single(seq1_rc, seq2_rc, n_0_error, n_1_error)
    )


def _overlap_graph(
    seqs: Sequence[str | Seq], n_0_error: int, n_1_error: int
) -> nx.Graph:
    """Create simple graph with seqs at vertices and edges when sequences overlap.

    Two sequences are considered overlapping if they share n_0_error bases with
    no errors or n_1_error bases with at most 1 error.

    Args:
        seqs: Sequence of sequences to build graph from
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch

    Returns:
        NetworkX Graph where vertex i corresponds to seqs[i]
    """
    seqs = [Seq(s) if isinstance(s, str) else s for s in seqs]

    if len(set(str(seq) for seq in seqs)) != len(seqs):
        warnings.warn("Input sequences to overlap_graph are not unique", UserWarning)

    g = nx.Graph()
    for i in range(len(seqs)):
        g.add_node(i)
    for i, seq1 in enumerate(seqs[:-1]):
        for j, seq2 in enumerate(seqs[i + 1 :], start=i + 1):
            if seqs_overlap(seq1, seq2, n_0_error, n_1_error):
                g.add_edge(i, j)

    return g


def overlap_inds(
    seqs: Sequence[str | Seq], subset_seqs: Dict[int, bool], n_0_error: int, n_1_error: int
) -> Dict[int, bool]:
    """Get indices of sequences that are in a connected component with any seed sequence.

    Computes indices of sequences in seqs that are in a connected component with any
    sequence containing any of the seed sequences (or its reverse complement) according
    to the overlap graph. Returns 0-based indices.

    Args:
        seqs: Sequences to analyze
        subset_seqs: Dict whose keys are indices into above list corresponding to sequences that
            have the seed sequence, and whose values are bools corresponding to whether the
            sequence is in forward orientation (True) or reverse orientation (False)
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch

    Returns:
        List of indices into seqs of sequences in connected components with any seed sequence
    """
    seqs = [Seq(s) if isinstance(s, str) else s for s in seqs]

    g = _overlap_graph(seqs, n_0_error, n_1_error)
    components: List[set] = nx.connected_components(g)

    # Find sequences containing seeds as substrings (is_subseq check)
    # Get the indices of all sequences that have a seed

    # Get all sets (connected components) of indices, where at least one index corresponds to a
    # seed-containing sequence
    seed_components: Dict[int, bool] = {}
    for component_set in components:
        for idx, is_reverse_complement in subset_seqs.items():
            if idx in component_set:
                seed_components.update({c: is_reverse_complement for c in component_set})
    return seed_components
