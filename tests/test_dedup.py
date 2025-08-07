"""Unit tests for outward_assembly.dedup module."""

import random

import networkx as nx
import pytest

from outward_assembly.dedup import (
    EMPTY_KMER_SENTINEL_HASH,
    ORIENT_STRICT,
    ORIENT_TOLERANT,
    MinimizerParams,
    ReadPair,
    _assign_to_buckets,
    _canonical_kmer,
    _extract_minimizer,
    _get_bucket_keys,
    _reverse_complement,
)


def _random_seq(length: int, rng: random.Random) -> str:
    """Generate random DNA sequence of specified length."""
    return "".join(rng.choices(["A", "C", "G", "T"], k=length))


class TestHelperFunctions:
    """Test sequence manipulation helper functions."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_reverse_complement_standard_bases(self):
        assert _reverse_complement("ACGT") == "ACGT"
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("TTTT") == "AAAA"
        assert _reverse_complement("GCGC") == "GCGC"

    @pytest.mark.fast
    @pytest.mark.unit
    def test_reverse_complement_with_n(self):
        assert _reverse_complement("ACGTN") == "NACGT"
        assert _reverse_complement("NNNNN") == "NNNNN"

    @pytest.mark.fast
    @pytest.mark.unit
    def test_reverse_complement_empty(self):
        assert _reverse_complement("") == ""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_canonical_kmer_lexicographic_selection(self):
        assert _canonical_kmer("AAAA") == "AAAA"  # AAAA vs TTTT
        assert _canonical_kmer("TTTT") == "AAAA"  # Same result
        assert _canonical_kmer("ACGT") == "ACGT"  # ACGT vs ACGT (palindrome)
        assert _canonical_kmer("AAAC") == "AAAC"  # AAAC vs GTTT
        assert _canonical_kmer("GTTT") == "AAAC"  # Same result


class TestMinimizerExtraction:
    """Test minimizer extraction functions."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_extract_minimizer_normal_window(self):
        params = MinimizerParams(num_windows=2, window_len=20, kmer_len=7)
        seq = "A" * 20 + "C" * 20  # 40bp sequence

        # Test first window (all A's - should give consistent result)
        hash1 = _extract_minimizer(seq, 0, params)
        hash2 = _extract_minimizer(seq, 0, params)
        assert hash1 == hash2  # Should be deterministic

        # Test second window (all C's - should give different result)
        hash3 = _extract_minimizer(seq, 1, params)
        assert hash1 != hash3  # Different sequences should give different hashes

    @pytest.mark.fast
    @pytest.mark.unit
    def test_extract_minimizer_with_n_bases(self):
        params = MinimizerParams(num_windows=1, window_len=10, kmer_len=3)
        seq_with_n = "AANAAAANAA"
        seq_without_n = "AAGAAAAGAA"

        hash_with_n = _extract_minimizer(seq_with_n, 0, params)
        hash_without_n = _extract_minimizer(seq_without_n, 0, params)

        # Should skip N-containing kmers and find valid ones
        assert hash_with_n != EMPTY_KMER_SENTINEL_HASH
        assert hash_without_n != EMPTY_KMER_SENTINEL_HASH

    @pytest.mark.fast
    @pytest.mark.unit
    def test_extract_minimizer_window_too_short(self):
        params = MinimizerParams(num_windows=1, window_len=10, kmer_len=7)
        seq = "AAAAA"  # 5bp sequence, need 7bp kmer

        hash_result = _extract_minimizer(seq, 0, params)
        assert hash_result == EMPTY_KMER_SENTINEL_HASH

    @pytest.mark.fast
    @pytest.mark.unit
    def test_extract_minimizer_sequence_too_short(self):
        "Collected windows longer than sequence, should succeed with a sentinel hash."
        params = MinimizerParams(num_windows=2, window_len=10, kmer_len=7)
        seq = "AAAAACCGGTT"  # 11bp sequence, second window is too short

        hash_result = _extract_minimizer(seq, 1, params)
        assert hash_result == EMPTY_KMER_SENTINEL_HASH

    @pytest.mark.fast
    @pytest.mark.unit
    def test_extract_minimizer_sequence_matches_window_matches_kmer(self):
        params = MinimizerParams(num_windows=1, window_len=11, kmer_len=11)
        seq = "AAAAACCGGTT"  # 11bp sequence

        hash_result = _extract_minimizer(seq, 0, params)
        assert hash_result != EMPTY_KMER_SENTINEL_HASH

    @pytest.mark.fast
    @pytest.mark.unit
    def test_extract_minimizer_all_N_window(self):
        params = MinimizerParams(num_windows=1, window_len=10, kmer_len=3)
        seq = "NNNNNNNNNN"

        hash_result = _extract_minimizer(seq, 0, params)
        assert hash_result == EMPTY_KMER_SENTINEL_HASH

    @pytest.mark.fast
    @pytest.mark.unit
    def test_get_bucket_keys(self):
        params = MinimizerParams(num_windows=2, window_len=20, kmer_len=7)
        rng = random.Random("hello")
        rp = ReadPair(
            "test", _random_seq(40, rng), _random_seq(40, rng), "I" * 40, "I" * 40
        )

        # Strict mode should generate num_windows² keys:
        # num_windows minimizers in the fwd seq x num_windows in the rev seq
        keys = _get_bucket_keys(rp, params, ORIENT_STRICT)
        assert len(keys) == 4
        assert all(isinstance(key, tuple) and len(key) == 2 for key in keys)

        # Tolerant mode should generate 2*num_windows² keys:
        # num_windows minimizers in the fwd seq x num_windows in the rev seq,
        # times 2 for swapping fwd/rev
        keys = _get_bucket_keys(rp, params, ORIENT_TOLERANT)
        assert len(keys) == 8
        assert all(isinstance(key, tuple) and len(key) == 2 for key in keys)


class TestBucketing:
    """Test bucketing functions."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_assign_to_buckets_correct_assignment(self):
        # Identical read pairs should go to same bucket
        params = MinimizerParams(num_windows=1, window_len=10, kmer_len=3)

        rp1 = ReadPair("read1", "A" * 10, "T" * 10, "I" * 10, "I" * 10)
        rp2 = ReadPair("read2", "A" * 10, "T" * 10, "I" * 10, "I" * 10)

        buckets = _assign_to_buckets([rp1, rp2], params, ORIENT_STRICT)

        # Should have at least one bucket containing both reads
        assert any(len(bucket_indices) >= 2 for bucket_indices in buckets.values())

    @pytest.mark.fast
    @pytest.mark.unit
    def test_assign_to_buckets_multiple_buckets(self):
        # Reads should appear in multiple buckets (since num_windows > 1)
        params = MinimizerParams(num_windows=2, window_len=10, kmer_len=3)
        rng = random.Random("hello")
        rp = ReadPair(
            "read1",
            _random_seq(20, rng),
            _random_seq(20, rng),
            "I" * 20,
            "I" * 20,
        )

        buckets = _assign_to_buckets([rp], params, ORIENT_STRICT)

        # Read should appear in multiple buckets (W² = 4 buckets)
        total_appearances = sum(0 in indices for indices in buckets.values())
        assert total_appearances >= 4

    @pytest.mark.fast
    @pytest.mark.unit
    def test_dups_share_buckets(self):
        """Test with realistic-ish data that two read pairs which are duplicates
        with a single base error in each of fwd/rev read share a bucket. This test
        is probabilistic, though with fixed random seed should be consistent."""
        rng = random.Random(1234)
        read_len = 150

        def _seq_with_error(seq):
            error_loc = rng.choice(range(read_len))
            if rng.random() < 0.5:  # sometimes offset
                start = 1
            else:
                start = 0
            return (
                seq[start:error_loc]
                + rng.choice(["A", "C", "G", "T"])
                + seq[error_loc + 1 :]
            )

        for rep in range(100):
            # create a pair of reads that are dups but not perfect dups
            rp1 = ReadPair(
                "pair1",
                _random_seq(read_len, rng),
                _random_seq(read_len, rng),
                "I" * read_len,
                "I" * read_len,
            )
            rp2_fwd = _seq_with_error(rp1.fwd_seq)
            rp2_rev = _seq_with_error(rp1.rev_seq)
            rp2 = ReadPair(
                "pair2", rp2_fwd, rp2_rev, "I" * len(rp2_fwd), "I" * len(rp2_rev)
            )
            buckets = _assign_to_buckets([rp1, rp2], MinimizerParams(), ORIENT_STRICT)
            # Should have at least one bucket containing both reads
            assert any(len(bucket_indices) == 2 for bucket_indices in buckets.values())


class TestReadPairClass:
    """Test ReadPair class functionality."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_post_init_sequences_uppercase(self):
        rp = ReadPair("test", "acgt", "tgca", "IIII", "IIII")
        assert rp.fwd_seq == "ACGT"
        assert rp.rev_seq == "TGCA"

    @pytest.mark.fast
    @pytest.mark.unit
    def test_mean_qual_calculation(self):
        # Phred 33: '!' = 0, 'I' = 40
        rp = ReadPair("test", "AAAA", "TTTT", "!!!!", "IIII")
        expected_mean = (0 + 0 + 0 + 0 + 40 + 40 + 40 + 40) / 8
        assert rp.mean_qual() == expected_mean

    @pytest.mark.fast
    @pytest.mark.unit
    def test_mean_qual_empty_qualities(self):
        rp = ReadPair("test", "AAAA", "TTTT", "", "")
        assert rp.mean_qual() == 0.0


class TestParameterValidation:
    """Test parameter validation."""

    @pytest.mark.fast
    @pytest.mark.unit
    def test_minimizer_params_kmer_too_large(self):
        with pytest.raises(ValueError, match="kmer_len .* must be <= window_len"):
            MinimizerParams(window_len=5, kmer_len=7)

    @pytest.mark.fast
    @pytest.mark.unit
    def test_minimizer_params_valid(self):
        params = MinimizerParams(window_len=10, kmer_len=5)
        assert params.window_len == 10
        assert params.kmer_len == 5
        # it's legal if surprising to have window and kmer length match
        params = MinimizerParams(window_len=47, kmer_len=47)
        assert params.window_len == 47
        assert params.kmer_len == 47
