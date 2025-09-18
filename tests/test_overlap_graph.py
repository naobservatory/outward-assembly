import networkx
import pytest

from outward_assembly.basic_seq_operations import SeqOrientation
from outward_assembly.overlap_graph import (
    _traverse_subgraph_and_orient,
    get_overlapping_sequence_ids,
)


@pytest.mark.fast
@pytest.mark.unit
def test_overlap_inds_multiple_seeds():
    """Test that overlap_inds correctly identifies sequences
    in connected components with seed-containing sequences."""
    # Set up test sequences
    seqs = [
        "ACGTACGTACGT",  # seq0: contains seed0
        "ACGTACGTTTTT",  # seq1: overlaps with seq0 but doesn't contain seed
        "GTTTTTCGTACGTG",  # seq2: overlaps with seq1 but not with seq0
        "GGGGGGTTTTCCCCTTT",  # seq3: contains rc(seed1)
        "AAAACCCCGGGAAAGG",  # seq4: overlaps with rc(seq3)
        "ZZZZZZZZZZZZ",  # seq5: no overlap, no seed (should be excluded)
    ]

    # Define seeds
    seq_ids_containing_seeds = {
        0: SeqOrientation.FORWARD,  # seq0: forward direction
        3: SeqOrientation.FORWARD,  # seq3: forward direction
    }

    # Test connected component behavior
    result = get_overlapping_sequence_ids(
        seqs, seq_ids_containing_seeds, n_0_error=5, n_1_error=7
    )

    # Should include sequences containing seeds and those connected to them
    assert sorted(result) == [0, 1, 2, 3, 4]


@pytest.mark.fast
@pytest.mark.unit
@pytest.mark.parametrize(
    "initial_orientation", (SeqOrientation.FORWARD, SeqOrientation.REVERSE)
)
def test_traverse_subgraph_and_orient(initial_orientation):
    """
    Tests that starting from an initial node, traverse_subgraph_and_orient finds all
    connected components and re-orients them all with respect to the initial node's
    orientation.

    The overlap graph looks like this, and we will orient all nodes with respect to node 2:
    0 --fwd-- 1 --rev-- 2 --fwd-- 3     5
                        |
                        --rev-- 4
    """
    # Construct overlap graph
    g = networkx.Graph()
    g.add_nodes_from([0, 1, 2, 3, 4, 5])
    g.add_edge(0, 1, orientation=SeqOrientation.FORWARD)
    g.add_edge(1, 2, orientation=SeqOrientation.REVERSE)
    g.add_edge(2, 3, orientation=SeqOrientation.FORWARD)
    g.add_edge(2, 4, orientation=SeqOrientation.REVERSE)

    oriented_connected_components = _traverse_subgraph_and_orient(
        g, reference_seq_idx=2, reference_seq_orientation=initial_orientation
    )
    expected_orientations = {
        0: SeqOrientation.REVERSE,
        1: SeqOrientation.REVERSE,
        2: SeqOrientation.FORWARD,
        3: SeqOrientation.FORWARD,
        4: SeqOrientation.REVERSE,
    }
    if initial_orientation == SeqOrientation.FORWARD:
        assert oriented_connected_components == expected_orientations
    else:
        assert oriented_connected_components == {
            node: -1 * orientation for node, orientation in expected_orientations.items()
        }
