#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=align_minimap2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --array=1-1000%20
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

mkdir -p logs

MINIMAP2="hpc quay.io/biocontainers/minimap2:2.28--h577a1d6_4 minimap2"
SAMTOOLS="hpc staphb/samtools:1.19 samtools"

# ---------------------------------------------------------
# Defaults for interactive execution
# ---------------------------------------------------------

THREADS="${SLURM_CPUS_PER_TASK:-4}"
TASK="${SLURM_ARRAY_TASK_ID:-1}"
OFFSET="${OFFSET:-0}"

TASK=$((TASK + OFFSET))

# ---------------------------------------------------------
# Convert task number into Set and Individual
# ---------------------------------------------------------

SET=$(( (TASK-1)/100 + 1 ))
IND=$(( (TASK-1)%100 + 1 ))

ROOT="/home/projects/btip26_divsim/artillum_test/third_test/generated_individuals_trial"

SETDIR="${ROOT}/set_${SET}"
READDIR="${SETDIR}/reads"
BAMDIR="${SETDIR}/bam"

mkdir -p "$BAMDIR"

# ---------------------------------------------------------
# Reference (same folder as this script)
# ---------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="/home/projects/btip26_divsim/artillum_test/third_test/00_reference.fasta"

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Cannot find reference:"
    echo "$REF"
    exit 1
fi

# ---------------------------------------------------------
# Build minimap2 index once
# ---------------------------------------------------------

if [[ ! -f "${REF}.mmi" ]]; then
    echo "Building minimap2 index..."
    $MINIMAP2 -d "${REF}.mmi" "$REF"
fi

# ---------------------------------------------------------
# Locate FASTQ files
# ---------------------------------------------------------

R1=$(find "$READDIR" -maxdepth 1 -name "${SET}_individual_${IND}_*_R1.fq")
R2=$(find "$READDIR" -maxdepth 1 -name "${SET}_individual_${IND}_*_R2.fq")

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: FASTQ files not found."
    echo "Set : $SET"
    echo "Individual : $IND"
    echo "Directory : $READDIR"
    exit 1
fi

PREFIX=$(basename "$R1" _R1.fq)

BAM="${BAMDIR}/${PREFIX}.bam"

echo "======================================"
echo "Task       : $TASK"
echo "Set        : $SET"
echo "Individual : $IND"
echo "Threads    : $THREADS"
echo "Reference  : $REF"
echo "R1         : $R1"
echo "R2         : $R2"
echo "Output BAM : $BAM"
echo "======================================"

# ---------------------------------------------------------
# Align and sort
# ---------------------------------------------------------

$MINIMAP2 \
    -t "$THREADS" \
    -ax sr \
    -R "@RG\tID:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA" \
    "$REF" \
    "$R1" \
    "$R2" \
| $SAMTOOLS sort \
    --threads "$THREADS" \
    -o "$BAM"

# ---------------------------------------------------------
# Index BAM
# ---------------------------------------------------------

$SAMTOOLS index \
    --threads "$THREADS" \
    "$BAM"

echo
echo "======================================"
echo "Finished!"
echo "Output:"
echo "$BAM"
echo "======================================"
