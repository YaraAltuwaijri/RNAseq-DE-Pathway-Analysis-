#!/bin/bash
#SBATCH --job-name=fastqc_job          # Job name
#SBATCH --output=../logs/slurm-%j.out  # Standard output and error log

# Project and directory variables
PROJECT_DIR="$(dirname "$(pwd)")"
INPUT_DIR=${1:-$PROJECT_DIR/data/raw}                 # First argument: Input directory (default: data/raw)
OUTPUT_DIR=${2:-$PROJECT_DIR/outputs/fastqc_results}  # Second argument: Output directory (default: outputs/fastqc_results)
THREADS=${3:-4}                                       # Third argument: Number of threads (default: 4)

# Load FastQC module (or activate its environment if using Conda)
module load fastqc/0.12.1-gcc-13.2.0

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC for all FASTQ files in the input directory
fastqc "$INPUT_DIR"/*.fastq -o "$OUTPUT_DIR" -t "$THREADS"

# Check exit status and notify
if [ $? -eq 0 ]; then
    echo "FastQC completed successfully for all files in $INPUT_DIR."
else
    echo "FastQC encountered an error. Check the logs for details."
fi
