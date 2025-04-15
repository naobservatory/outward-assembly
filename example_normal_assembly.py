from outward_assembly import s3_files_with_prefix, outward_assembly

def main():
    # --- User-defined parameters ---
    # Directory where the results will be stored locally
    data_dir = "<path to store results locally>"
    # Path to the seed sequence file
    seed_path = "<path to seed sequence file>"
    # Name of the S3 bucket where the SIZ files are stored
    bucket_name = "<bucket name>"
    # List of S3 prefixes for the SIZ files (e.g., ["folder/file_prefix_1", "folder/file_prefix_2"])
    prefixes = ["<path(s) to siz file(s)>"]
    # --- End of user-defined parameters ---

    # Define the working directory
    workdir = data_dir + "/work"

    # Get all S3 paths
    paths = [
        p for prefix in prefixes for p in s3_files_with_prefix(bucket_name, prefix)
    ]

    _ = outward_assembly(
        s3_paths=paths,
        seed_path=seed_path,
        output_path=workdir + "/final_output.fasta",
        adapters_path="./default_adapters.fasta",
        work_dir_parent=workdir,
    )

if __name__ == "__main__":
    main()

