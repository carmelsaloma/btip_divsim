#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=chr1_chr2_divcalc
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

mkdir -p logs

CHROM_DIRECT="/home/projects/btip26_divsim/chromosomes_final"
OUTPUT_FILE="/home/projects/btip26_divsim/diversity_summary_chr1_chr2.tsv"
SCRIPT_DIRECT="/home/projects/btip26_divsim"

[[ -d "$CHROM_DIRECT" ]] || { echo "ERROR: chromosome directory not found: $CHROM_DIRECT" >&2; exit 1; }

if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found in PATH." >&2
    exit 1
fi

if ! python3 -c "import Bio, numpy" &> /dev/null; then
    echo "ERROR: biopython and/or numpy not importable. Activate the right conda env first (e.g. 'conda activate <env>')." >&2
    exit 1
fi

echo "Calculating nucleotide diversity (Nei & Li 1979) for chr1 and chr2 only in $CHROM_DIRECT"

python3 "$SCRIPT_DIRECT/calculate_diversity.py" \
  --chrom-dir "$CHROM_DIRECT" \
  --pattern 'chr[12].fasta' \
  --out "$OUTPUT_FILE"

echo "Done. Diversity summary (chr1 + chr2 only) written to: $OUTPUT_FILE"
