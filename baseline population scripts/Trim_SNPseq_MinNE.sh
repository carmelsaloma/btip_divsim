#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=trimming_nonneutral
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=1-200%10
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --partition=batch

MIN_NE=25000
NE_FILE="ne_values.txt"

# Let SLURM handle the loop tracking dynamically via its Array Task ID
i=$SLURM_ARRAY_TASK_ID

HALF_NE=$(sed -n "${i}p" "$NE_FILE")
NE=$((2 * HALF_NE))

FILE="runs/Ne_${NE}_task${i}/run_${NE}/run_${NE}_SNPseq.csv"
OUTFILE="runs/Ne_${NE}_task${i}/run_${NE}/trimmed_SNPseq_${NE}.csv"

if [[ -f "$FILE" ]]; then
    {
        head -n 1 "$FILE"
        tail -n +2 "$FILE" | head -n "$MIN_NE"
    } > "$OUTFILE"
    echo "Task $i - Created: $OUTFILE"
else
    echo "Task $i - Missing: $FILE"
    exit 1
fi

