#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=generate_fasta
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=100-200%10
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --partition=batch

# Ensure log directory exists
mkdir -p logs

NE_LIST=ne_values.txt

# Extract the half-Ne value matching the current Slurm task ID
HALF_NE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$NE_LIST")

echo "$HALF_NE done processing."

NE=$((2 * HALF_NE))

echo "$NE done processing."

# Define the run directory using the exact nested structure requested
RUNDIR="runs/Ne_${NE}_task${SLURM_ARRAY_TASK_ID}/run_${NE}"

# Define input and output file paths relative to the root directory
indices_file="${RUNDIR}/run_${NE}_indices.txt"
ref_file="${RUNDIR}/run_${NE}_tempseq.fa"
out_file="run_${NE}.fasta"

# Dynamic fallback to handle either naming convention for the SNP CSV
if [[ -f "${RUNDIR}/trimmed_SNPseq_${NE}.csv" ]]; then
    csv_file="${RUNDIR}/trimmed_SNPseq_${NE}.csv"
else
    csv_file="${RUNDIR}/run_${NE}_SNPseq.csv"
fi

echo "Processing Array Task ${SLURM_ARRAY_TASK_ID} | NE = ${NE}"

# Execute the python script with the required arguments

source ~/miniconda3/etc/profile.d/conda.sh
conda activate py3env

python3 generate_fasta2.py \
    -i "${indices_file}" \
    -r "${ref_file}" \
    -s "${csv_file}" \
    -o "${out_file}" \
    --one-based
