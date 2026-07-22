#!/usr/bin/env bash

inputfasta="merged_nonneutral.fasta"
outputdir="per_sample"

mkdir -p "$outputdir"

gawk -v outputdir="$outputdir" '
  BEGIN {
    sample_count = 0
  }
  /^>/ {
    # Extract sample ID (between >1_ and chr)
    if (match($0, /^>1_([0-9]+)chr/, arr)) {
      sample_id = arr[1]

      # Zero-pad to 5 digits
      file_id = sprintf("%05d", sample_id)

      # Build output filename
      file = outputdir "/nonneutral_block_sample_" file_id ".fasta"

      # If this is the first header for this sample, increment counter and echo progress
      if (!(sample_id in seen)) {
        sample_count++
        seen[sample_id] = 1
        printf("Processing sample %s -> %s\n", sample_id, file) > "/dev/stderr"
      }
    } else {
      printf("WARNING: header did not match expected pattern, skipping: %s\n", $0) > "/dev/stderr"
      file = ""
    }
  }
  {
    if (file != "") {
      print > file
    }
  }
  END {
    printf("\nSplit complete: %d samples written to folder %s\n", sample_count, outputdir) > "/dev/stderr"
  }
' "$inputfasta"


