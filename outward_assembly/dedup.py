"""
Post-assembly read-pair deduplication with error-tolerant matching.

This module provides deduplication of read pairs using minimizer-based
bucket for efficiency. Deduplication is tolerant to small alignment
shifts and sequencing errors.
"""

import sys
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from zlib import crc32
from typing import Dict, List, Literal, Optional, Set, Tuple

import networkx as nx

# No Bio import; we represent sequences as strings rather than Bio.Seq:
# - Don't need any Bio machinery
# - Python string operations are faster than the corresponding Seq operations

ORIENT_STRICT = "strict"
ORIENT_TOLERANT = "tolerant"


@dataclass
class MinimizerParams:
    """Minimizer configuration (rarely needs changing)."""

    num_windows: int = 2  # Number of windows per read
    window_len: int = 20  # Base pairs per window
    kmer_len: int = 7  # K-mer size for minimizers

    def __post_init__(self):
        if self.kmer_len > self.window_len:
            raise ValueError(
                f"kmer_len ({self.kmer_len}) must be <= window_len ({self.window_len})"
            )


@dataclass
class ReadPair:
    """Container for a read pair with deduplication support."""

    read_id: str
    fwd_seq: str
    rev_seq: str
    fwd_qual: str
    rev_qual: str
    exemplar_id: Optional[str] = field(default=None, init=False)

    def __post_init__(self):
        """Ensure sequences are uppercase."""
        self.fwd_seq = self.fwd_seq.upper()
        self.rev_seq = self.rev_seq.upper()

    def mean_qual(self) -> float:
        """Calculate mean Phred quality across both reads."""
        quals = [ord(c) - 33 for c in self.fwd_qual + self.rev_qual]
        return sum(quals) / len(quals) if quals else 0.0


##
# Assign read pairs to buckets based on minimizers.
# Each read pair will be assigned to multiple buckets.
# With high probability, duplicate read pairs will be assigned to at least
# one bucket in common, so we only need to do all-against-all read pair comparisons
# within each bucket.
##

# Complement table for reverse complement (including N)
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def _canonical_kmer(kmer: str) -> str:
    """Return lexicographically smaller of kmer and its reverse complement."""
    rc = _reverse_complement(kmer)
    return min(kmer, rc)

def _hash_kmer(kmer: str) -> int:
    """Hash a kmer to an int. The actual hash used is an implementation detail,
    but the result must be stable run-to-run (so no default Python hash)."""
    return crc32(kmer.encode())

def _extract_minimizer(seq: str, window_idx: int, params: MinimizerParams) -> int:
    """
    Extract the minimizer hash from a specific window of the sequence.

    Args:
        seq: DNA sequence
        window_idx: Which window to process (0-based)
        params: Minimizer parameters

    Returns:
        Hash of the lexicographically smallest canonical k-mer in the window
    """
    start = window_idx * params.window_len
    end = min(len(seq), start + params.window_len)

    if end - start < params.kmer_len:
        # Window too short to contain a k-mer - return consistent hash
        return _hash_kmer("EMPTY")

    # Find minimizer (smallest hash) in this window
    bigger_than_hash = sys.maxsize + 1
    min_hash = bigger_than_hash
    for i in range(start, end - params.kmer_len + 1):
        kmer = seq[i : i + params.kmer_len]
        if "N" not in kmer:  # Skip k-mers with ambiguous bases
            canonical = _canonical_kmer(kmer)
            h = _hash_kmer(canonical)
            if h < min_hash:
                min_hash = h

    return min_hash if min_hash != bigger_than_hash else _hash_kmer("EMPTY")


def _get_bucket_keys(
    read_pair: ReadPair, params: MinimizerParams, orientation: str
) -> Set[Tuple[int, int]]:
    """
    Generate all bucket keys for a read pair based on minimizers.

    Returns set of (forward_hash, reverse_hash) tuples that serve as bucket keys.
    """
    # Extract minimizers from each window
    fwd_hashes = [
        _extract_minimizer(read_pair.fwd_seq, i, params)
        for i in range(params.num_windows)
    ]
    rev_hashes = [
        _extract_minimizer(read_pair.rev_seq, i, params)
        for i in range(params.num_windows)
    ]

    # Generate all hash pairs
    keys = {(fh, rh) for fh in fwd_hashes for rh in rev_hashes}

    # In tolerant mode, also consider swapped orientation
    if orientation == ORIENT_TOLERANT:
        keys |= {(rh, fh) for fh in fwd_hashes for rh in rev_hashes}

    return keys


def _assign_to_buckets(
    read_pairs: List[ReadPair], minimizer_params: MinimizerParams, orientation: str
) -> Dict[Tuple[int, int], List[int]]:
    """Assign read pairs to buckets by minimizers. Return a Dict {bucket_key : indices}
    where bucket_key is a tuple of ints (kmer hashes) and indices are relative to the
    input list of read pairs."""
    buckets = defaultdict(list)
    for idx, rp in enumerate(read_pairs):
        keys = _get_bucket_keys(rp, minimizer_params, orientation)
        for key in keys:
            buckets[key].append(idx)
    return buckets

