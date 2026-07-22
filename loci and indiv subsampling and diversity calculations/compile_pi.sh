#!/bin/bash
#SBATCH --account=btip26_divsim
#SBATCH --job-name=compile_pi_indivandloci
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# 06_compile_pi_indivandloci.sh
#
# Compiles all pixy_results/set_${SET}/pixy_set${SET}_loci${LOCI}_indiv${N}_pi.txt
# files produced by Script 05 into a single combined text file.

mkdir -p logs

# ---------------------------
# Config
# ---------------------------
SET_LIST=($(seq 1 50))
LOCI_LIST=(30 60 90 120 150 180 210 240 270 300)
INDIV_LIST=(1 2 3 4 5 10 20 30 50 100)

IN_DIRECT="pixy_results_final"
OUT_DIRECT="pixy_compiled_final"
OUT_FILE="${OUT_DIRECT}/final_compiled_pi.txt"

mkdir -p "$OUT_DIRECT"

# ---------------------------
# Guards
# ---------------------------
[[ -d "$IN_DIRECT" ]] || { echo "ERROR: $IN_DIRECT not found. Run Script 05 first." >&2; exit 1; }

# ---------------------------
# Compile
# ---------------------------
first=true
n_compiled=0

echo "=== Aggregating pixy results across 50 sets ==="

for SET in "${SET_LIST[@]}"; do

    SET_DIR="${IN_DIRECT}/set_${SET}"

    if [[ ! -d "$SET_DIR" ]]; then
        echo "WARNING: $SET_DIR not found. Skipping set ${SET}." >&2
        continue
    fi

    for LOCI in "${LOCI_LIST[@]}"; do
        for N in "${INDIV_LIST[@]}"; do

            pi_file="${SET_DIR}/pixy_set${SET}_loci${LOCI}_ind${N}_pi.txt"

            if [[ ! -f "$pi_file" ]]; then
                echo "WARNING: $pi_file not found. Skipping." >&2
                continue
            fi

            if $first; then
                # write the combined header once, with extra leading columns
                awk -v set="$SET" -v loci="$LOCI" -v n="$N" \
                    'NR==1 { print "set\tloci\tn_indiv\t" $0; next } { print set"\t"loci"\t"n"\t"$0 }' \
                    "$pi_file" > "$OUT_FILE"
                first=false
            else
                # subsequent files: skip header line, append data rows
                awk -v set="$SET" -v loci="$LOCI" -v n="$N" \
                    'NR==1 { next } { print set"\t"loci"\t"n"\t"$0 }' \
                    "$pi_file" >> "$OUT_FILE"
            fi

            n_compiled=$((n_compiled + 1))

        done
    done
    echo "  compiled set ${SET}"
done

if $first; then
    echo "ERROR: No pixy pi result files were found under $IN_DIRECT -- nothing to compile." >&2
    exit 1
fi

echo "Done. Compiled ${n_compiled} files. Final output written to: $OUT_FILE"
