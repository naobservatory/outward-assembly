from outward_assembly import s3_files_with_prefix, outward_assembly
from outward_assembly.kmer_freq_filter import _high_freq_kmers_split_files
from outward_assembly.io_helpers import process_s3_paths
import os
from dotenv import load_dotenv

load_dotenv()

def main():
    # Define your S3 prefixes
    data_dir = "<path to store results locally>"
    workdir = data_dir + "/work"

    # Prefixes of paths to reads on s3 for siz files, for example, <folder>/<file prefix>
    prefixes = ["<path to siz file>"]

    # Get all S3 paths
    paths = [
        p for prefix in prefixes for p in s3_files_with_prefix("<bucket name>", prefix)
    ]

    # OPTIONAL: Run high frequency k-mer filter if the file doesn't exist
    if not os.path.exists(workdir + "/kmers/high_freq_kmers.fasta"):
        high_freq_path = _high_freq_kmers_split_files(
            process_s3_paths([paths[i] for i in range(1, len(paths), 100)]),
            workdir,
            min_kmer_freq=1000,
            num_parallel=10,
        )
    else:
        high_freq_path = workdir + "/kmers/high_freq_kmers.fasta"

    _ = outward_assembly(
        s3_paths=paths,
        seed_path=data_dir + "/seed_seq.fasta",
        warm_start_path=data_dir + "/warm_start.fasta",
        output_path=workdir + "/final_output.fasta",
        adapters_path="./default_adapters.fasta",
        work_dir_parent=workdir,
        high_freq_kmers_path=high_freq_path,
        read_subset_k=25,
        use_batch=True,
        batch_queue="<batch queue name>",
        batch_workdir="<s3 working path>",
    )

if __name__ == "__main__":
    main()

