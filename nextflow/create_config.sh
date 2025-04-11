#!/bin/bash
# Usage:
#   ./create_config.sh [s3_dir] [sample_path_file] [ref_fasta_path] [bbduk_k] [batch-queue] [tower_token (optional)]
#
# Example:
#   ./create_config.sh "s3://my-bucket/path" "/my/project/data/sample_paths.txt" "/my/project/data/query_kmers.fasta" 27 "my-queue" "my-tower-token"

# Set default values if not provided on the command line
S3_DIR=${1}
SAMPLE_PATH_FILE=${2}
REF_FASTA_PATH=${3}
BBDUK_K=${4}
BATCH_QUEUE=${5}
TOWER_TOKEN=${6:-} 

if [ -z "${S3_DIR}" ] || [ -z "${SAMPLE_PATH_FILE}" ] || [ -z "${REF_FASTA_PATH}" ] || [ -z "${BBDUK_K}" ] || [ -z "${BATCH_QUEUE}" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: ./create_config.sh [s3_dir] [sample_path_file] [ref_fasta_path] [bbduk_k] [batch-queue] [tower_token (optional)]"
    exit 1
fi

# Determine the directory of the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cat <<EOF > "${SCRIPT_DIR}/nextflow.config"
params {
  // You can override the base_dir here
  base_dir = "${S3_DIR}"
  s3_files = "${SAMPLE_PATH_FILE}"
  ref_fasta_path = "${REF_FASTA_PATH}"

  // bbduk params
  kmer = ${BBDUK_K}
}

includeConfig "\${projectDir}/config/output.config"
includeConfig "\${projectDir}/config/logging.config"
includeConfig "\${projectDir}/config/resources.config"
includeConfig "\${projectDir}/config/profiles.config"
includeConfig "\${projectDir}/config/containers.config"
process.queue = '${BATCH_QUEUE}'
tower.accessToken = '${TOWER_TOKEN}'
EOF
