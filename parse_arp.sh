#!/bin/bash

# This is the parser script for the .arp file to get the loci & actual SNP sequences per sample. 

# Check if input is ok

echo
echo "

=========================================

        popgenART: ARP parsing 

=========================================
"
echo 

while getopts "p:a:" opt; do
  case "$opt" in

    p) prefix="$OPTARG" ;;
    a) arpfile="$OPTARG" ;;

    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [[ -z "$prefix" || -z "$arpfile" ]]; then
  echo "Usage: $0 -p <prefix> -a <arpfile>" >&2
  exit 1
fi

echo "Complete arguments for arlequin parsing."





indices_array=${prefix}_indices.txt
echo "" > $indices_array

echo "Saving indices representing the loci for the SNPs."

mapfile -t line_nums < <(grep -n "polymorphic positions on chromosome" "$arpfile" | cut -d: -f1)

chrom_num=1
for line_num in "${line_nums[@]}"; do
	echo "Chrom ---- ${chrom_num}" >> "$indices_array"
	target_line=$((line_num + 1))
	line=$(sed -n "${target_line}p" "$arpfile")
	clean_line=$(echo "$line" | tr -d '#' | xargs)
	IFS=',' read -ra polyInd <<< "$clean_line"

	for i in "${polyInd[@]}"; do
    		val=$(echo "$i" | xargs)
    		echo "$val" >> "$indices_array"
  	done
	
	((chrom_num++))
	
done


echo "Indices are saved in ${indices_array}"

## Get the SNP sequence per sample and save to a csv file.

echo "Saving SNP sequences in the csv file."

awk '
  # Start flag when inside SampleData block
  /SampleData= *\{/ {inBlock=1; next}

  # End flag when closing bracket of SampleData block
  /^\}/ && inBlock {inBlock=0; next}

  # When inside SampleData, print lines with 3 fields (Sample lines)
  inBlock && NF==3 {print $1","$2","$3}
' ${arpfile} > ${prefix}_SNPseq.csv

(echo "SampleID,Col2,Sequence"; awk '
  /SampleData= *\{/ {inBlock=1; next}
  /^\}/ && inBlock {inBlock=0; next}
  inBlock && NF==3 {print $1","$2","$3}
' ${arpfile}) > ${prefix}_SNPseq.csv

echo "Saved the SNPs to a csv file."


