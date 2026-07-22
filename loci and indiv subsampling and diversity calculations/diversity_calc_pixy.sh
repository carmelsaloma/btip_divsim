#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=calc_pi_pixy
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --reservation=btip26
#SBATCH --array=3-50
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -eo pipefail

mkdir -p logs

# ---------------------------
# Environment
# ---------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pixy-env

# ---------------------------
# Config
# ---------------------------
PIXY=pixy
LOCUS_SIZE=25000  # 25 kb per locus
THREADS=${SLURM_CPUS_PER_TASK:-4}

LOCI_LIST=(30 60 90 120 150 180 210 240 270 300)
INDIV_LIST=(1 2 3 4 5 10 20 30 50 100)

# Current array task = simulation set
SET=${SLURM_ARRAY_TASK_ID}

# Directories
VCF_DIR="/home/projects/btip26_divsim/artillum_test/third_test/final_subsampled_individuals"
OUT_DIR="pixy_results_final/set_${SET}"

mkdir -p "$OUT_DIR"

# Tool for extracting sample names if pop file is missing
BCFTOOLS="hpc quay.io/biocontainers/bcftools:1.22--h3a4d415_1 bcftools"

echo "=========================================="
echo "Measuring pi for SET ${SET}"
echo "=========================================="

for LOCI in "${LOCI_LIST[@]}"; do

    # Calculate window size dynamically based on the locus count
    WINDOW_SIZE=$(( LOCI * LOCUS_SIZE ))

    for N in "${INDIV_LIST[@]}"; do

        tag="set${SET}_loci${LOCI}_ind${N}"
        vcf_file="${VCF_DIR}/${tag}.vcf.gz"

        if [[ ! -f "$vcf_file" ]]; then
            echo "WARNING: $vcf_file not found. Skipping."
            continue
        fi

	# -------------------------------------------------------------
        # Guard: Check for VCF Index (.tbi or .csi) -- Required by pixy
        # -------------------------------------------------------------
        if [[ ! -f "${vcf_file}.tbi" && ! -f "${vcf_file}.csi" ]]; then
            echo "Indexing VCF: $vcf_file"
            $BCFTOOLS index -t "$vcf_file"
        fi

        # Population file path inside OUT_DIR
        pop_file="${OUT_DIR}/pop_indiv${N}.txt"

        # Generate population file
        $BCFTOOLS query -l "$vcf_file" | awk '{print $1"\tpop1"}' > "$pop_file"       

        out_prefix="pixy_${tag}"
        expected_out="${OUT_DIR}/${out_prefix}_pi.txt"

        if [[ -f "$expected_out" ]]; then
            echo "Output $expected_out already exists. Skipping."
            continue
        fi

        echo "=== Measuring pi: Set ${SET} | Loci ${LOCI} | Indiv ${N} ==="

        $PIXY --stats pi \
            --vcf "$vcf_file" \
            --populations "$pop_file" \
            --window_size "$WINDOW_SIZE" \
            --n_cores "$THREADS" \
            --output_folder "$OUT_DIR" \
            --output_prefix "$out_prefix"

        echo "  wrote: $expected_out"

    done
done

echo "Done. All pixy results for SET ${SET} are written to: $OUT_DIR"
