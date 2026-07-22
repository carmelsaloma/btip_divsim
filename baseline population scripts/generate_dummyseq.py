#!/usr/bin/env python3

"""
popgenART: Generate a dummy reference FASTA from an fsc2-style parameter CSV.

Fixes two bugs in the bash version:
  - No out-of-bounds array indexing (which silently produced empty-string
    "bases" and shortened sequences below the requested length).
  - Base pool is not accidentally shared/accumulated across chromosomes.

Reads Num_Chrom, Num_Loci_1a, and GC_Con from the CSV and writes exactly
Num_Loci_1a bases per chromosome, GC-weighted, straight to FASTA.

Usage:
  python3 generate_ref.py -c params.csv -p prefix
"""
import argparse
import csv
import random
import sys


def read_params(path):
    params = {}
    with open(path) as f:
        reader = csv.reader(f)
        next(reader, None)  # header: Parameter,Value
        for row in reader:
            if not row or row[0].startswith("<EndOfParFile>"):
                continue
            if len(row) >= 2 and row[0].strip():
                params[row[0].strip()] = row[1].strip()
    return params


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--csv", required=True, help="fsc2 parameter CSV")
    ap.add_argument("-p", "--prefix", required=True, help="output prefix")
    ap.add_argument("--seed", type=int, default=None, help="random seed (optional, for reproducibility)")
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    params = read_params(args.csv)

    try:
        num_chrom = int(params["Num_Chrom"])
        num_loci = int(params["Num_Loci_1a"])
        gc_con = float(params["GC_Con"])
    except KeyError as e:
        sys.exit(f"ERROR: missing required parameter in CSV: {e}")

    if not (0.0 <= gc_con <= 1.0):
        sys.exit(f"ERROR: GC_Con must be between 0 and 1, got {gc_con}")

    gc = gc_con / 2.0       # proportion each for G and C
    at = (1.0 - gc_con) / 2.0  # proportion each for A and T
    bases = ["A", "T", "G", "C"]
    weights = [at, at, gc, gc]

    print(f"Number of chromosomes: {num_chrom}", file=sys.stderr)
    print(f"Loci per chromosome: {num_loci}", file=sys.stderr)
    print(f"GC content: {gc_con} (A={at:.3f} T={at:.3f} G={gc:.3f} C={gc:.3f})", file=sys.stderr)

    out_path = f"{args.prefix}_tempseq.fa"
    with open(out_path, "w") as out:
        for chrom in range(1, num_chrom + 1):
            seq = "".join(random.choices(bases, weights=weights, k=num_loci))
            out.write(f">Chromosome {chrom}\n{seq}\n")
            if chrom % 100 == 0 or chrom == num_chrom:
                print(f"  generated {chrom}/{num_chrom}", file=sys.stderr)

    print(f"Done: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()$


