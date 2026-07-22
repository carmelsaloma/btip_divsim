#!/bin/bash
#SBATCH --job-name=subset_breadth
#SBATCH --partition=batch
#SBATCH --account=btip26_divsim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-50%10
#SBATCH --output=subset_breadth_%A_%a.out
#SBATCH --error=subset_breadth_%A_%a.err

set -euo pipefail

BCFTOOLS="quay.io/biocontainers/bcftools:1.22--h3a4d415_1"

# Current array task = simulation set
SET=${SLURM_ARRAY_TASK_ID}


# Genome structure
LOCUS_SIZE=25000
TOTAL_LOCI=300


# Breadth levels (10% increments)
BREADTHS=(
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


echo
echo "========================================"
echo "Processing set ${SET}"
echo "========================================"



SETDIR="generated_individuals_trial/set_${SET}"

VCFDIR="${SETDIR}/vcf"

OUTDIR="${SETDIR}/final_subsampling"


mkdir -p "${OUTDIR}"


INPUT="${VCFDIR}/all100.vcf.gz"



# Generate fixed locus coordinates

LOCUS_FILE="${SETDIR}/loci_coordinates.txt"

> "${LOCUS_FILE}"


for LOCUS in $(seq 1 ${TOTAL_LOCI}); do

    START=$(( (LOCUS - 1) * LOCUS_SIZE + 1 ))
    END=$(( LOCUS * LOCUS_SIZE ))

    echo -e "locus${LOCUS}\t${START}\t${END}" >> "${LOCUS_FILE}"

done



for NLOCI in "${BREADTHS[@]}"; do


    echo
    echo "Creating ${NLOCI} loci subset"


    REGION_FILE="${OUTDIR}/regions_loci${NLOCI}.txt"


    awk -v n="${NLOCI}" '
    NR <= n {
        print "1:"$2"-"$3
    }
    ' "${LOCUS_FILE}" > "${REGION_FILE}"



    OUTPUT="${OUTDIR}/set${SET}_loci${NLOCI}.vcf.gz"



    hpc ${BCFTOOLS} \
        bcftools view \
        -R "${REGION_FILE}" \
        -Oz \
        -o "${OUTPUT}" \
        "${INPUT}"



    hpc ${BCFTOOLS} \
        bcftools index \
        -f "${OUTPUT}"



done



rm -f "${LOCUS_FILE}"


echo
echo "========================================"
echo "Completed set ${SET}"
echo "========================================"
