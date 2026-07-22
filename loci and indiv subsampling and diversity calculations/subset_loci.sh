#!/bin/bash

set -euo pipefail


# ==========================================
# Final loci breadth subsampling
# SCREEN VERSION
#
# Sets:
#   1-50
#
# Breadth:
#   30,60,90,120,150,180,210,240,270,300 loci
#
# Sequential genomic blocks
#
# Each locus:
#   25,000 bp
#
# ==========================================


BCFTOOLS="quay.io/biocontainers/bcftools:1.22--h3a4d415_1"


# -------------------------------
# Settings
# -------------------------------

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


WINDOW=25000


CHR="1_9170_chr1"



# -------------------------------
# Main loop
# -------------------------------

for SET in "${SETS[@]}"
do


    echo
    echo "====================================="
    echo "Processing SET ${SET}"
    echo "====================================="



    SETDIR="generated_individuals_trial/set_${SET}"


    INPUT="${SETDIR}/vcf/all100.vcf.gz"


    OUTDIR="${SETDIR}/final_subsampling"



    if [ ! -f "${INPUT}" ]; then

        echo "WARNING: Missing input"
        echo "${INPUT}"
        echo "Skipping SET ${SET}"

        continue

    fi



    mkdir -p "${OUTDIR}"



    for LOCI in "${LOCI_LEVELS[@]}"
    do


        REGION_FILE="${OUTDIR}/regions_loci${LOCI}.txt"


        OUTPUT="${OUTDIR}/set${SET}_loci${LOCI}.vcf.gz"



        END=$((LOCI * WINDOW))



        echo
        echo "SET ${SET}"
        echo "Creating loci ${LOCI}"
        echo "Region:"
        echo "${CHR}:1-${END}"



        # Generate sequential loci regions

        awk \
        -v chr="${CHR}" \
        -v loci="${LOCI}" \
        -v size="${WINDOW}" '

        BEGIN{

            for(i=0;i<loci;i++){

                start=(i*size)+1

                end=(i+1)*size

                print chr "\t" start "\t" end

            }

        }' > "${REGION_FILE}"



        echo "Extracting variants..."



        hpc ${BCFTOOLS} bcftools view \
            -R "${REGION_FILE}" \
            -Oz \
            -o "${OUTPUT}" \
            "${INPUT}"



        echo "Indexing..."



        hpc ${BCFTOOLS} bcftools index \
            -f \
            "${OUTPUT}"



        rm -f "${REGION_FILE}"



        echo "Finished SET ${SET} loci ${LOCI}"


    done


done



echo
echo "====================================="
echo "ALL 50 SETS COMPLETED"
echo "====================================="
