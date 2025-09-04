import pytest

from outward_assembly.overlap_graph import overlap_inds


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
        0: True,  # seq0: forward direction
        3: True,  # seq3: forward direction
    }

    # Test connected component behavior
    result = overlap_inds(seqs, seq_ids_containing_seeds, n_0_error=5, n_1_error=7)

    # Should include sequences containing seeds and those connected to them
    assert sorted(result) == [0, 1, 2, 3, 4]
