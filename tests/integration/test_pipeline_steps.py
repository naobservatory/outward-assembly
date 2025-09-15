import os
import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from outward_assembly.pipeline_steps import _subset_contigs


@pytest.fixture
def temp_workdir_with_contigs():
    """Create a temporary directory with mock contigs for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)

        # Create megahit output directory
        megahit_dir = workdir / "megahit_out_iter1-1"
        megahit_dir.mkdir(parents=True)

        # Create mock contigs
        contigs = [
            SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="contig1"),
            SeqRecord(Seq("TTTTTTTTTTTTTTTTT"), id="contig2"),
            SeqRecord(Seq("GGGGGGAAAAAACCCCCC"), id="contig3"),
            SeqRecord(Seq("TACGCGTACGCGTACGTCGTA"), id="contig4"),
            SeqRecord(Seq("ATAGACGATGGCGAGGCGCGTA"), id="contig5"),
        ]

        # Write contigs to file
        contigs_path = megahit_dir / "final.contigs.fa"
        SeqIO.write(contigs, contigs_path, "fasta")

        yield workdir


@pytest.mark.fast
@pytest.mark.unit
def test_subset_contigs_multiple_seeds(temp_workdir_with_contigs):
    """Test _subset_contigs with multiple seed sequences."""
    workdir = temp_workdir_with_contigs

    # Define multiple seeds
    seed_seqs = [
        Seq("ACGTAC"),  # seed1: in contig1, RC in contig4
        Seq("AAAAAA"),  # seed2: in contig3, RC in contig2
    ]

    # Run subset_contigs
    _subset_contigs(workdir, iter=1, seed_seqs=seed_seqs)

    # Check output
    filtered_path = workdir / "megahit_out_iter1-1" / "contigs_filtered.fasta"
    filtered_contigs = list(SeqIO.parse(filtered_path, "fasta"))

    # Should find 4 contigs (1, 2, 3, and 4)
    assert len(filtered_contigs) == 4

    # Check contig IDs
    contig_ids = [rec.id for rec in filtered_contigs]
    assert "contig1" in contig_ids
    assert "contig2" in contig_ids
    assert "contig3" in contig_ids
    assert "contig4" in contig_ids
    assert "contig5" not in contig_ids


@pytest.mark.fast
@pytest.mark.unit
def test_subset_contigs_with_overlaps(temp_workdir_with_contigs):
    """Test _subset_contigs with include_overlaps parameter."""
    workdir = temp_workdir_with_contigs

    # Define seeds
    seed_seqs = [Seq("ATATAAAGGCC")]

    # Create an additional contig that overlaps with seed but doesn't contain it
    megahit_dir = workdir / "megahit_out_iter1-1"
    contigs_path = megahit_dir / "final.contigs.fa"
    # Add a contig with overlap to existing contigs
    with open(contigs_path, "a") as f:
        f.write(">contig6\nGCGCGTAATATAAAGGCC\n")  # Contains seed and overlaps with contig5

    # Run subset_contigs with include_overlaps=True
    _subset_contigs(workdir, iter=1, seed_seqs=seed_seqs, include_overlaps=True)

    # Check output
    filtered_path = megahit_dir / "contigs_filtered.fasta"
    filtered_contigs = list(SeqIO.parse(filtered_path, "fasta"))

    # Should include contig5 due to overlap
    contig_ids = [rec.id for rec in filtered_contigs]
    assert "contig5" in contig_ids
