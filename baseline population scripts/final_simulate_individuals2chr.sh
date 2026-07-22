#!/bin/bash

# --- CONFIGURATION ---
SAMTOOLS="hpc staphb/samtools:1.19 samtools"
set -euo pipefail

INPUT_DIRECT="/home/projects/btip26_divsim/chromosomes_final"
OUTPUT_DIRECT="/home/projects/btip26_divsim/generated_individuals_trial2"
POPULATION_FILE="${OUTPUT_DIRECT}/population.txt"
IND="$SLURM_ARRAY_TASK_ID"

mkdir -p "$OUTPUT_DIRECT"

# ---------- 1. Generate Population File (Only once) ----------
# This block ensures all jobs share the same 100 IDs.
# If the file exists, it will NOT be overwritten (preserving consistency).
if [[ ! -f "$POPULATION_FILE" ]]; then
    exec 200>"${OUTPUT_DIRECT}/.pop.lock"
    flock -x 200
    if [[ ! -f "$POPULATION_FILE" ]]; then
        echo "Generating unique population list (100 IDs)..."
        shuf -i 1-10000 -n 100 > "$POPULATION_FILE"
    fi
    exec 200>&-
fi

# ---------- 2. Get the Sample ID for this Individual ----------
# Pulls the N-th line from the file where N is the SLURM task ID
SAMPLE_ID=$(sed -n "${IND}p" "$POPULATION_FILE")

# ---------- 3. Build the Individual (Overwrites existing file) ----------
OUTFILE="${OUTPUT_DIRECT}/individual_${IND}.fasta"

# The '>' operator clears the file if it exists, ensuring a fresh, clean write.
> "$OUTFILE"

# Loop through all 10 chromosomes, using the SAME Sample ID for all
for chr in {1..2}; do
    SEQID="1_${SAMPLE_ID}_chr${chr}"
    
    # Append the sequence to the clean output file
    $SAMTOOLS faidx "${INPUT_DIRECT}/chr${chr}.fasta" "$SEQID" >> "$OUTFILE"
done

echo "Individual ${IND} built successfully using Sample ID: ${SAMPLE_ID}"
