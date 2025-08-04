#!/bin/bash
#SBATCH --job-name=download_and_convert_job   # Job name
#SBATCH --output=../logs/slurm-%j.out  # Standard output

# Set variables from arguments
PROJECT_DIR="$(dirname "$(pwd)")"             # Path to the main project directory
OPTION_FILE=$1                                # First argument: option file (e.g., SRR_Acc_list.txt)
OUTPUT_DIR=${2:-$PROJECT_DIR/data/raw}        # Second argument: output directory (default: 'data/raw')
MEMORY=${3:-8G}                               # Third argument: memory (e.g., 8G)
THREADS=${4:-4}                               # Fourth argument: number of threads (default: 4)

#Loading SRA-tools module for running fasterq. Can also be done using conda
module load sra-tools/3.0.3-gcc-13.2.0

# Create the output directory if it doesn't already exist
mkdir -p "$OUTPUT_DIR"

# Loop through each accession in the option file
while IFS= read -r accession; do
    echo "Downloading and converting $accession to FASTQ..."

    # Run fasterq-dump
    fasterq-dump --split-files "$accession" -O "$OUTPUT_DIR" --threads "$THREADS"

    # Check if the download and conversion were successful
    if [ $? -eq 0 ]; then
        echo "$accession downloaded and converted successfully."
    else
        echo "Failed to download and convert $accession."
    fi
done < "$OPTION_FILE"
