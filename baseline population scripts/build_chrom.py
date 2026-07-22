#!/usr/bin/env python3
"""
Build chromosome-level FASTA files by stitching together per-locus sequences
for each sample, using a map.txt file that lists (in order) which loci belong
to which chromosome.

Uses Biopython (Bio.SeqIO) to read and write FASTA files instead of manual
parsing.

Directory assumptions:
    neutral block per-sample fastas:     /home/projects/btip26_divsim/neutral_block/per_sample
    non-neutral block per-sample fastas: /home/projects/btip26_divsim/nonneutral_block/per_sample

Each per-sample fasta is expected to be named:
    neutral_block_{sample}.fasta
    nonneutral_block_{sample}.fasta

Each fasta header inside those files is expected to look like:
    >1_{sample}_chr{loci}

map.txt is expected to look like this (a summary line followed by one line of
comma-separated loci for that chromosome), repeated for each chromosome:

    chr1: total=150  neutral=139  nonneutral=11
    10, 15, 67, 70, 75, 76, ...

    chr2: total=150  neutral=139  nonneutral=11
    4, 11, 18, 19, 25, ...

The order of loci within each comma-separated line = the order they get
concatenated in to build that chromosome's sequence.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_map(map_path):
    """Read map.txt -> {chromosome: [loci_id, loci_id, ...]} preserving list order.

    Expects blocks like:
        chr1: total=150  neutral=139  nonneutral=11
        10, 15, 67, 70, ...
    """
    chrom_loci = defaultdict(list)
    header_re = re.compile(r"^(chr\d+):")

    with open(map_path) as f:
        current_chrom = None
        for line_num, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line:
                continue

            match = header_re.match(line)
            if match:
                # this is a "chrN: total=... neutral=... nonneutral=..." line
                current_chrom = match.group(1)
                continue

            if current_chrom is None:
                # stray line before any chromosome header (e.g. title lines) - skip
                continue

            # this should be the comma-separated loci line for current_chrom
            loci = [x.strip() for x in line.split(",") if x.strip()]
            if not loci:
                print(f"WARNING: map.txt line {line_num} had no loci parsed, skipping", file=sys.stderr)
                continue
            chrom_loci[current_chrom].extend(loci)
            current_chrom = None  # reset until next header line is seen

    return chrom_loci


def load_fasta(path):
    """Parse a fasta file into {header: sequence} using Biopython's SeqIO.

    Note: SeqIO.parse strips the leading '>' and only keeps the part of the
    header up to the first whitespace as record.id. Since our headers
    (>1_{sample}_{loci}) have no spaces, record.id is the full label minus '>'.
    We add the '>' back here so lookups elsewhere can stay consistent.
    """
    seqs = {}
    for record in SeqIO.parse(path, "fasta"):
        seqs[f">{record.id}"] = str(record.seq)
    return seqs


def build_sample_dict(sample, neutral_dir, nonneutral_dir):
    """Load + merge the neutral and non-neutral fasta for one sample into one dict."""
    padded = f"{int(sample):05d}"
    neutral_path = os.path.join(neutral_dir, f"neutral_block_sample_{padded}.fasta")
    nonneutral_path = os.path.join(nonneutral_dir, f"nonneutral_block_sample_{padded}.fasta")

    sample_dict = {}

    if os.path.isfile(neutral_path):
        sample_dict.update(load_fasta(neutral_path))
    else:
        print(f"WARNING: missing neutral file for sample {sample}: {neutral_path}", file=sys.stderr)

    if os.path.isfile(nonneutral_path):
        sample_dict.update(load_fasta(nonneutral_path))
    else:
        print(f"WARNING: missing non-neutral file for sample {sample}: {nonneutral_path}", file=sys.stderr)

    return sample_dict


def build_chromosome_sequence(sample, chrom, loci_list, sample_dict):
    """Concatenate the sequence for each locus (in map.txt order) to build one chromosome."""
    pieces = []
    missing = 0
    for loci in loci_list:
        key = f">1_{sample}chr{loci}"
        seq = sample_dict.get(key)
        if seq is None:
            missing += 1
            continue
        pieces.append(seq)
    if missing:
        print(f"  WARNING: sample {sample} {chrom}: {missing} loci missing from sample_dict", file=sys.stderr)
    return "".join(pieces)


def parse_sample_list(samples_arg):
    """Accepts a range like '1-25000' or a comma-separated list like '1,2,3'."""
    samples = []
    for part in samples_arg.split(","):
        part = part.strip()
        if "-" in part:
            start, end = part.split("-")
            samples.extend(str(s) for s in range(int(start), int(end) + 1))
        else:
            samples.append(part)
    return samples


def main():
    parser = argparse.ArgumentParser(description="Build chromosome-level fasta files per sample (Biopython version).")
    parser.add_argument("--map", default="map.txt", help="Path to map.txt (chromosome, loci columns)")
    parser.add_argument("--neutral-dir", default="/home/projects/btip26_divsim/neutral_block/per_sample")
    parser.add_argument("--nonneutral-dir", default="/home/projects/btip26_divsim/nonneutral_block/per_sample")
    parser.add_argument("--out-dir", default="./chromosomes_final")
    parser.add_argument("--samples", required=True, help="e.g. '1-25000' or '1,2,3,4'")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("Parsing map.txt ...")
    chrom_loci = parse_map(args.map)
    chromosomes = sorted(chrom_loci.keys())
    print(f"Found {len(chromosomes)} chromosomes in map.txt: {chromosomes}")

    samples = parse_sample_list(args.samples)
    print(f"Processing {len(samples)} samples ...")

    # collect SeqRecords per chromosome, then write each chromosome file once at the end
    chrom_records = {chrom: [] for chrom in chromosomes}

    for i, sample in enumerate(samples, start=1):
        sample_dict = build_sample_dict(sample, args.neutral_dir, args.nonneutral_dir)

        for chrom in chromosomes:
            loci_list = chrom_loci[chrom]
            final_seq = build_chromosome_sequence(sample, chrom, loci_list, sample_dict)
            record = SeqRecord(Seq(final_seq), id=f"1_{sample}_{chrom}", description="")
            chrom_records[chrom].append(record)

        if i % 100 == 0 or i == len(samples):
            print(f"  processed {i}/{len(samples)} samples")

    for chrom in chromosomes:
        out_path = os.path.join(args.out_dir, f"{chrom}.fasta")
        SeqIO.write(chrom_records[chrom], out_path, "fasta")
        print(f"  wrote {len(chrom_records[chrom])} records -> {out_path}")

    print(f"Done. Chromosome fastas written to: {args.out_dir}")


if __name__ == "__main__":
    main()
