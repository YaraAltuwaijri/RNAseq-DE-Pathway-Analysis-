#!/bin/bash
#SBATCH --output=../logs/slurm-%j.out  # Standard output and error log

# Define base project directory
baseDirectory="$(dirname "$(pwd)")" # Root directory
genomeDirectory="$baseDirectory/data/raw/ncbi_dataset/data/GCF_000146045.2"

# Input and output directories

# Directory with input files (default to $PROJECT_DIR/data/processed)
inputDirectory=${1:-$baseDirectory/data/processed}
# Directory to store the alignment output (default to $PROJECT_DIR/outputs)
outputDirectory=${2:-$baseDirectory/outputs}
# Default memory (8G if not specified)
MEMORY=${3:-8G}
 # Default number of threads (4 if not specified)
THREADS=${4:-4}

# creating ouput directory if it doesn't exist
mkdir -p "$outputDirectory"

# STAR, samtools, and featureCounts containers
starContainer="$baseDirectory/containers/star-alignment_latest.sif"
samtoolsContainer="$baseDirectory/containers/samtools_latest.sif"
subreadContainer="$baseDirectory/containers/subread_latest.sif"

# Looping through all FASTQ files in input directory
for inputFile in "$inputDirectory"/*_trimmed.fastq; do

    fileName=$(basename "$inputFile")
    sampleID="${fileName%%_*}"

    # Defining results directory for this sample
    resultsDirectory="$outputDirectory/alignment_$sampleID"
    mkdir -p "$resultsDirectory"

    echo "Processing sample: $sampleID"
    echo "Results directory: $resultsDirectory"

    # STAR alignment and generating the BAM file
    echo "Running STAR alignment for $sampleID..."
    singularity exec --bind "$resultsDirectory","$baseDirectory" "$starContainer" \
        STAR --runThreadN "$THREADS" \
             --genomeDir "$genomeDirectory" \
             --readFilesCommand cat \
             --readFilesIn "$inputFile" \
             --outFileNamePrefix "$resultsDirectory/theAlignment" \
             --outSAMtype BAM SortedByCoordinate

    # Indexing the BAM file with samtools
    echo "Indexing BAM file for $sampleID..."
    singularity exec --bind "$resultsDirectory","$baseDirectory" "$samtoolsContainer" \
        samtools index "$resultsDirectory/theAlignmentAligned.sortedByCoord.out.bam"

    # Performing feature counting with subread
    echo "Counting features for $sampleID..."
    singularity exec --bind "$resultsDirectory","$baseDirectory" "$subreadContainer" \
        featureCounts -p -g gene_id -a "$genomeDirectory/genomic.gtf" \
                      -o "$resultsDirectory/counts.txt" "$resultsDirectory/theAlignmentAligned.sortedByCoord.out.bam"

    echo "Completed processing for sample: $sampleID"
    echo "----------------------------------------"
done

echo "All samples have been processed successfully."
