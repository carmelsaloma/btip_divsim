#!/bin/bash

set -euo pipefail

# ==========================================================
# 08_final_subset_individuals_screen.sh
#
# Purpose:
#   Nested individual subsampling
#
# Screen version
#
# Input:
#   generated_individuals_trial/set_X/final_subsampling/
#
# Output:
#   generated_individuals_trial/set_X/final_individual_subsampling/
#
# Sets:
#   1-50
#
# Loci:
#   30,60,90,120,150,180,210,240,270,300
#
# Individuals:
#   1,2,3,4,5,10,20,30,50,100
#
# ==========================================================


BCFTOOLS="quay.io/biocontainers/bcftools:1.22--h3a4d415_1"


SETS=( {1..50} )

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
echo "=========================================================="
echo "08 FINAL INDIVIDUAL SUBSAMPLING"
echo "=========================================================="
echo


for SET in "${SETS[@]}"
do

    echo
    echo "##########################################################"
    echo "SET ${SET}"
    echo "##########################################################"

    SETDIR="generated_individuals_trial/set_${SET}"

    INPUTDIR="${SETDIR}/final_subsampling"

    OUTDIR="${SETDIR}/final_individual_subsampling"

    LISTDIR="${OUTDIR}/sample_lists"

    mkdir -p "${OUTDIR}"
    mkdir -p "${LISTDIR}"


    REFERENCE_VCF="${INPUTDIR}/set${SET}_loci30.vcf.gz"

    if [ ! -f "${REFERENCE_VCF}" ]; then

        echo
        echo "WARNING"
        echo "Missing:"
        echo "${REFERENCE_VCF}"
        echo "Skipping SET ${SET}"
        echo

        continue

    fi


    echo
    echo "Generating sample lists..."


    ALL_SAMPLES="${LISTDIR}/all_samples.txt"

    hpc ${BCFTOOLS} \
        bcftools query \
        -l \
        "${REFERENCE_VCF}" \
        > "${ALL_SAMPLES}"


    TOTAL=$(wc -l < "${ALL_SAMPLES}")

    if [ "${TOTAL}" -ne 100 ]; then

        echo
        echo "ERROR"
        echo "Expected 100 samples."
        echo "Found ${TOTAL}"
        exit 1

    fi


    for N in "${IND_LEVELS[@]}"
    do

        head -n "${N}" "${ALL_SAMPLES}" \
            > "${LISTDIR}/ind${N}.txt"

    done


    echo "Sample lists complete."



    for LOCI in "${LOCI_LEVELS[@]}"
    do

        echo
        echo "----------------------------------------------------------"
        echo "SET ${SET} | LOCI ${LOCI}"
        echo "----------------------------------------------------------"

        INPUTVCF="${INPUTDIR}/set${SET}_loci${LOCI}.vcf.gz"

        if [ ! -f "${INPUTVCF}" ]; then

            echo "Missing:"
            echo "${INPUTVCF}"
            echo "Skipping."

            continue

        fi

        for N in "${IND_LEVELS[@]}"
        do

            SAMPLEFILE="${LISTDIR}/ind${N}.txt"

            OUTPUTVCF="${OUTDIR}/set${SET}_loci${LOCI}_ind${N}.vcf.gz"

            echo
            echo "Individuals: ${N}"

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

            echo "Finished:"
            echo "$(basename "${OUTPUTVCF}")"

        done

        echo
        echo "Completed loci ${LOCI}"

    done

    echo
    echo "##########################################################"
    echo "SET ${SET} COMPLETE"
    echo "##########################################################"

done


echo
echo "=========================================================="
echo "ALL 50 SETS FINISHED"
echo "=========================================================="
