#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=varcall_bcftools
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --array=1-50%20
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

mkdir -p logs

BCFTOOLS="hpc quay.io/biocontainers/bcftools:1.22--h3a4d415_1 bcftools"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
SET="${SLURM_ARRAY_TASK_ID:-1}"

ROOT="/home/projects/btip26_divsim/artillum_test/third_test/generated_individuals_trial"

SETDIR="${ROOT}/set_${SET}"
BAMDIR="${SETDIR}/bam"
VCFDIR="${SETDIR}/vcf"

mkdir -p "$VCFDIR"

ROOT="/home/projects/btip26_divsim/artillum_test/third_test"

REF="${ROOT}/00_reference.fasta"

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference not found:"
    echo "$REF"
    exit 1
fi

BAMLIST="${BAMDIR}/bamlist.txt"

find "$BAMDIR" -maxdepth 1 -name "*.bam" | sort > "$BAMLIST"

NUM_BAMS=$(wc -l < "$BAMLIST")

if [[ "$NUM_BAMS" -ne 100 ]]; then
    echo "ERROR: Expected 100 BAM files."
    echo "Found: $NUM_BAMS"
    exit 1
fi

echo "======================================"
echo "Set        : $SET"
echo "Reference  : $REF"
echo "BAM files  : $NUM_BAMS"
echo "Output     : ${VCFDIR}/all100.vcf.gz"
echo "======================================"

$BCFTOOLS mpileup \
    --threads "$THREADS" \
    -Ou \
    -f "$REF" \
    -b "$BAMLIST" \
|
$BCFTOOLS call \
    --threads "$THREADS" \
    -m \
    -A \
    -Oz \
    -o "${VCFDIR}/all100.vcf.gz"

$BCFTOOLS index \
    "${VCFDIR}/all100.vcf.gz"

echo
echo "======================================"
echo "Finished Set ${SET}"
echo "VCF:"
echo "${VCFDIR}/all100.vcf.gz"
echo "======================================"
