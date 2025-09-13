import warnings
from collections import deque
from typing import Dict, Sequence

import networkx as nx
from Bio.Seq import Seq

from outward_assembly.basic_seq_operations import SeqOrientation


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
    seq1: str | Seq,
    seq2: str | Seq,
    n_0_error: int,
    n_1_error: int,
    relative_orientation: SeqOrientation,
) -> bool:
    """
    Return whether two sequences overlap by at least n_0_error bases with no errors
    or n_1_error bases with at most 1 error. Detects reverse complement overlaps.

    Args:
        seq1: First sequence to check for overlap
        seq2: Second sequence to check for overlap
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch
        relative_orientation: The relative orientations of the sequences to check. If
            Forward, checks whether the sequences overlap in forward orientation relative
            to each other. If Reverse, checks whether the sequences overlap in the reverse
            compliment relative to each other.

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

    if relative_orientation is SeqOrientation.FORWARD:
        return _seqs_overlap_single(
            seq1, seq2, n_0_error, n_1_error
        ) or _seqs_overlap_single(seq1_rc, seq2_rc, n_0_error, n_1_error)
    else:
        return _seqs_overlap_single(
            seq1, seq2_rc, n_0_error, n_1_error
        ) or _seqs_overlap_single(seq1_rc, seq2, n_0_error, n_1_error)


def _overlap_graph(seqs: Sequence[str | Seq], n_0_error: int, n_1_error: int) -> nx.Graph:
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
            # Annotate overlap orientation as an attribute on the edge
            if seqs_overlap(seq1, seq2, n_0_error, n_1_error, SeqOrientation.FORWARD):
                g.add_edge(i, j, orientation=SeqOrientation.FORWARD)
            elif seqs_overlap(seq1, seq2, n_0_error, n_1_error, SeqOrientation.REVERSE):
                g.add_edge(i, j, orientation=SeqOrientation.REVERSE)
    return g


def _traverse_subgraph_and_orient(
    g: nx.Graph, seed_contig_idx: int, seed_orientation: SeqOrientation
) -> Dict[int, SeqOrientation]:
    """
    Given a graph where the nodes are contigs and edges represent connected (overlapping)
    contigs, find all nodes in a connected component with seed_contig_idx, and return
    them all in the proper orientation with respect to seed_contig_idx. The seed contig's
    neighbors are traversed via BFS.

    Args:
        g: A networkx graph whose nodes represent contigs. For every pair of contigs A and
           B, there is an edge between A and B iff the sequences of A and B overlap, and
           the edge is annotated with the metadata key "orientation", whose value is the
           relative orientation of A and B (forward or reverse)
        seed_contig_idx: The contig to return all connected components for
        seed_orientation: The initial orientation of the seed contig
    """
    orientations: Dict[int, SeqOrientation] = {seed_contig_idx: seed_orientation}
    # Queue of nodes whose neighbors we want to traverse next in the breadth-first search
    bfs_queue = deque([seed_contig_idx])
    visited = {seed_contig_idx}

    while bfs_queue:
        current = bfs_queue.popleft()
        current_orientation = orientations[current]
        for neighbor in g.neighbors(current):
            if neighbor not in visited:
                # The relative orientation between current and neighbor is stored as an
                # attribute on the edge between them
                neighbor_relative_orientation: SeqOrientation = g[current][neighbor][
                    "orientation"
                ]
                neighbor_canonical_orientation = (
                    current_orientation * neighbor_relative_orientation
                )
                orientations[neighbor] = neighbor_canonical_orientation
                visited.add(neighbor)
                bfs_queue.append(neighbor)

    return orientations


def overlap_inds(
    seqs: Sequence[str | Seq],
    subset_seqs: Dict[int, SeqOrientation],
    n_0_error: int,
    n_1_error: int,
) -> Dict[int, SeqOrientation]:
    """Get indices of sequences that are in a connected component with any seed sequence.

    Computes indices of sequences in seqs that are in a connected component with any
    sequence containing any of the seed sequences (or its reverse complement) according
    to the overlap graph. Returns 0-based indices.

    Args:
        seqs: Sequences to analyze
        subset_seqs: Dict to use to subset the above the sequences. The keys of this dict are
            indices of seqs, corresponding to sequences that contain a seed, and the values
            are the orientation of each sequence relative to the seed it contains
        n_0_error: Minimum overlap length for exact match
        n_1_error: Minimum overlap length when allowing 1 mismatch

    Returns:
        List of indices into seqs of sequences in connected components with any seed sequence
    """
    seqs = [Seq(s) if isinstance(s, str) else s for s in seqs]

    # Construct a graph where nodes are contigs and edges represent overlaps between
    # contigs
    g = _overlap_graph(seqs, n_0_error, n_1_error)

    all_connected_contigs: Dict[int, SeqOrientation] = {}

    # Search for contigs that are connected to seed-containing contigs
    for seed_contig_idx, seed_contig_orientation in subset_seqs.items():
        if seed_contig_idx in all_connected_contigs:
            # We've already seen this contig before in a previous connected component, so
            # don't process it again
            continue

        # Get all contigs that are connected to this seed-containing contig, along with
        # their orientation relative to this seed-containing contig
        oriented_connected_contigs = _traverse_subgraph_and_orient(
            g, seed_contig_idx=seed_contig_idx, seed_orientation=seed_contig_orientation
        )
        all_connected_contigs.update(oriented_connected_contigs)

    return all_connected_contigs
