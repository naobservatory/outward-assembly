import time
from pathlib import Path

import pytest

from outward_assembly.pipeline import _compute_assembly_metrics


@pytest.mark.fast
@pytest.mark.unit
def test_compute_assembly_metrics_empty_iterations():
    """Test that _compute_assembly_metrics handles empty inner_iterations correctly."""
    workdir = Path("/test/dir")
    inner_iterations = []
    start_time = time.time() - 100  # 100 seconds ago

    # Call the function with empty inner_iterations
    metrics = _compute_assembly_metrics(workdir, inner_iterations, start_time)

    # Check that all required fields are present with default values
    assert metrics["total_time"] >= 100
    assert metrics["inner_iterations"] == []
    assert metrics["final_contig_count"] == 0
    assert metrics["final_longest_contig_length"] == 0
    assert metrics["final_total_contig_length"] == 0
    assert metrics["final_read_pair_count"] == 0
    assert metrics["work_dir"] == str(workdir)
