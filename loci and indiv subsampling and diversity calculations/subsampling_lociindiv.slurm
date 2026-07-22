#!/bin/bash
#SBATCH --job-name=08_final_subset_individuals
#SBATCH --partition=batch
#SBATCH --account=btip26_divsim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --reservation=btip26
#SBATCH --array=3-50%20
#SBATCH --output=logs/08_subset_individuals_%A_%a.out
#SBATCH --error=logs/08_subset_individuals_%A_%a.err

set -euo pipefail

BCFTOOLS="quay.io/biocontainers/bcftools:1.22--h3a4d415_1"

SET=${SLURM_ARRAY_TASK_ID}

SETDIR="generated_individuals_trial/set_${SET}"

INPUTDIR="${SETDIR}/final_subsampling"

OUTDIR="${SETDIR}/final_individual_subsampling"

LISTDIR="${OUTDIR}/sample_lists"

mkdir -p "${OUTDIR}"
mkdir -p "${LISTDIR}"

LOCI_LEVELS=(
30
60
90
120
150
180
210
240
270
300
)

IND_LEVELS=(
1
2
3
4
5
10
20
30
50
100
)

echo
echo "========================================"
echo "SET ${SET}"
echo "========================================"

REFERENCE_VCF="${INPUTDIR}/set${SET}_loci30.vcf.gz"

if [ ! -f "${REFERENCE_VCF}" ]; then
    echo "ERROR: Missing ${REFERENCE_VCF}"
    exit 1
fi

ALL_SAMPLES="${LISTDIR}/all_samples.txt"

if [ ! -f "${ALL_SAMPLES}" ]; then

    hpc ${BCFTOOLS} \
        bcftools query \
        -l \
        "${REFERENCE_VCF}" \
        > "${ALL_SAMPLES}"

fi

TOTAL=$(wc -l < "${ALL_SAMPLES}")

if [ "${TOTAL}" -ne 100 ]; then
    echo "ERROR: Expected 100 samples, found ${TOTAL}"
    exit 1
fi

for N in "${IND_LEVELS[@]}"
do

    SAMPLEFILE="${LISTDIR}/ind${N}.txt"

    if [ ! -f "${SAMPLEFILE}" ]; then

        head -n "${N}" "${ALL_SAMPLES}" \
            > "${SAMPLEFILE}"

    fi

done

for LOCI in "${LOCI_LEVELS[@]}"
do

    INPUTVCF="${INPUTDIR}/set${SET}_loci${LOCI}.vcf.gz"

    if [ ! -f "${INPUTVCF}" ]; then

        echo "Missing ${INPUTVCF}"

        continue

    fi

    echo
    echo "LOCI ${LOCI}"

    for N in "${IND_LEVELS[@]}"
    do

        SAMPLEFILE="${LISTDIR}/ind${N}.txt"

        OUTPUTVCF="${OUTDIR}/set${SET}_loci${LOCI}_ind${N}.vcf.gz"

        if [ -f "${OUTPUTVCF}" ]; then

            echo "Skipping $(basename "${OUTPUTVCF}")"

            continue

        fi

        echo "Creating $(basename "${OUTPUTVCF}")"

        hpc ${BCFTOOLS} \
            bcftools view \
            -S "${SAMPLEFILE}" \
            -Oz \
            -o "${OUTPUTVCF}" \
            "${INPUTVCF}"

        hpc ${BCFTOOLS} \
            bcftools index \
            -f \
            "${OUTPUTVCF}"

    done

done

echo
echo "SET ${SET} COMPLETE"
