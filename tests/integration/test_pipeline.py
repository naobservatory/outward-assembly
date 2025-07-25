import logging
import os
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from dotenv import load_dotenv

from outward_assembly.pipeline import outward_assembly

load_dotenv()


@pytest.mark.slow
@pytest.mark.e2e
@pytest.mark.requires_tools
@pytest.mark.requires_aws
@pytest.mark.parametrize(
    "use_batch,seed_config",
    [
        (True, {"file": "seed_seq.fasta", "name": "single_seed"}),
        (False, {"file": "seed_seq.fasta", "name": "single_seed"}),
        (False, {"file": "multi_seed_seq.fasta", "name": "multi_seed"}),
    ],
    ids=["batch-single-seed", "local-single-seed", "local-multi-seed"],
)
def test_pipeline(temp_workdir, use_batch, seed_config):
    """End to end test of full pipeline, following a simulated genome with
    structure ABCBD -- note the repeated B -- with three combinations:
    1. Batch mode with single seed
    2. Local mode with single seed
    3. Local mode with multiple seeds
    """
    # Capture real time output since it's a long test
    handler = logging.StreamHandler(sys.stdout)
    logging.getLogger().addHandler(handler)
    logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Extract seed file and scenario name from the seed_config
        seed_file = seed_config["file"]
        scenario_name = seed_config["name"]

        logging.info(
            f"Starting pipeline test: {scenario_name} in {'batch' if use_batch else 'local'} mode"
        )

        s3_paths = [
            "s3://nao-testing/outward-assembly-test-data/pipeline-end-to-end/simulated-abcbd-reads.fastq.zst"
        ]
        data_dir = Path(__file__).parent.parent / "data/pipeline-end-to-end/"

        # Define the output path
        output_path = temp_workdir / f"final_output_{scenario_name}.fasta"

        kwargs = {
            "s3_paths": s3_paths,
            "seed_path": data_dir / seed_file,
            "output_path": output_path,
            "max_iters": 2,
            "high_freq_kmers_path": data_dir / "fake_hfk.fasta",
            "freq_filter_k": 23,
            "work_dir_parent": temp_workdir,
            "adapters_path": data_dir / "fake_adapters.fasta",
            "warm_start_path": data_dir / "fake_warm_start.fasta",
        }

        if use_batch:
            kwargs.update(
                {
                    "use_batch": True,
                    "batch_workdir": os.getenv("BATCH_WORKDIR"),
                    "batch_queue": os.getenv("BATCH_QUEUE"),
                    "tower_token": os.getenv("TOWER_ACCESS_TOKEN"),
                }
            )

        outward_assembly(**kwargs)

        # Verify output exists
        assert output_path.exists(), f"Expected output file {output_path} not found"

        # Verify contigs were generated
        contigs = list(SeqIO.parse(output_path, "fasta"))
        assert len(contigs) > 0, "No contigs were generated"

        logging.info(
            f"Finished pipeline test: {scenario_name} in {'batch' if use_batch else 'local'} mode"
        )
    finally:
        logging.getLogger().removeHandler(handler)
