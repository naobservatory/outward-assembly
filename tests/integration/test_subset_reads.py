import os
from pathlib import Path

import pytest
from dotenv import load_dotenv

from outward_assembly.io_helpers import _count_lines, process_s3_paths
from outward_assembly.pipeline_steps import (
    _subset_split_files_batch,
    _subset_split_files_local,
)

load_dotenv()


# Runs the test twice, once with "use_batch=True" and once with "use_batch=False", where use_batch indicates whether to use batch mode
@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.requires_tools
@pytest.mark.requires_aws
@pytest.mark.parametrize("use_batch", [True, False])
def test_subset_reads(temp_workdir, use_batch):
    """Integration test of _subset_split_files.
    We pass in two split files with two read pairs each. The dummy contigs match the first
    forward read and last reverse read, so we should get two read pairs out.
    Tests both local and batch processing paths.
    """
    s3_records = process_s3_paths(
        [
            "s3://nao-testing/outward-assembly-test-data/read-filtering/test_reads_div0.fastq.zst",
            "s3://nao-testing/outward-assembly-test-data/read-filtering/test_reads_div1.fastq.zst",
        ]
    )
    data_dir = Path(__file__).parent.parent / "data/read-filtering/"

    kwargs = {
        "s3_records": s3_records,
        "ref_fasta_path": data_dir / "fake_contigs.fasta",
        "workdir": temp_workdir,
        "read_subset_k": 21,
    }

    if use_batch:
        kwargs.update(
            {
                "batch_workdir": os.getenv("BATCH_WORKDIR"),
                "batch_queue": os.getenv("BATCH_QUEUE"),
                "tower_token": os.getenv("TOWER_ACCESS_TOKEN"),
            }
        )
        _subset_split_files_batch(**kwargs)
    else:
        _subset_split_files_local(**kwargs)

    assert _count_lines(temp_workdir / "reads_1.fastq") == 8
    assert _count_lines(temp_workdir / "reads_2.fastq") == 8
