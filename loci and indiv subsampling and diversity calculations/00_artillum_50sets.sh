#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=artillum_50sets
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --array=1-1000%20
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

mkdir -p logs

INPUT_ROOT="/home/projects/btip26_divsim/artillum_test/third_test/generated_individuals_trial"

# ---------------------------------------------------------
# Convert array job (1-5000) into:
#   Set         = 1-50
#   Individual  = 1-100
# ---------------------------------------------------------

OFFSET=${OFFSET:-0}
TASK=$((SLURM_ARRAY_TASK_ID + OFFSET))

SET=$(( (TASK - 1) / 100 + 1 ))
IND=$(( (TASK - 1) % 100 + 1 ))

INPUT_DIR="${INPUT_ROOT}/set_${SET}"
OUTPUT_DIR="${INPUT_DIR}/reads"

mkdir -p "$OUTPUT_DIR"

# Verify that the set directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Directory not found:"
    echo "  $INPUT_DIR"
    exit 1
fi

# Locate the correct FASTA (there should only be one)
INPUT_FASTA=$(find "$INPUT_DIR" \
    -maxdepth 1 \
    -name "${SET}_individual_${IND}_*.fasta" \
    -print -quit)

if [[ -z "$INPUT_FASTA" ]]; then
    echo "ERROR: Could not find FASTA for:"
    echo "  Set ${SET}"
    echo "  Individual ${IND}"
    exit 1
fi

PREFIX=$(basename "$INPUT_FASTA" .fasta)

echo "=========================================="
echo "Task        : ${TASK}"
echo "Set         : ${SET}"
echo "Individual  : ${IND}"
echo "Input FASTA : ${INPUT_FASTA}"
echo "=========================================="

hpc arboradmin/art-illumina:250429 art_illumina \
    -na \
    -i "$INPUT_FASTA" \
    -l 150 \
    -f 10 \
    -nf 1 \
    -p \
    -m 200 \
    -s 10 \
    -ir 0.0001 \
    -ir2 0.00015 \
    -dr 0.00011 \
    -dr2 0.00022 \
    -qs 0 \
    -qs2 0 \
    -ss HS25 \
    -o "${OUTPUT_DIR}/${PREFIX}_R"

echo "Finished ${PREFIX}"
