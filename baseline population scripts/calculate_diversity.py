#!/usr/bin/env python3
"""
Calculate nucleotide diversity (pi) per Nei & Li (1979):

    pi = sum_{i<j} k_ij / (n(n-1)/2)

where k_ij is the number of nucleotide differences between sequences i and j,
and n is the number of sequences in the sample.

Rather than doing all n(n-1)/2 pairwise sequence comparisons directly (which
does not scale -- e.g. 10,000 sequences = ~50 million pairs per chromosome),
this uses the standard site-based reformulation: at each site, if a
nucleotide occurs with count c_a (a in {A,C,G,T}) out of n_valid sequences
with a called base at that site, the number of *differing* pairs at that
site is:

    n_valid*(n_valid-1)/2 - sum_a c_a*(c_a-1)/2

Summing this across all sites gives sum_{i<j} k_ij directly, in O(n*L) time
instead of O(n^2*L).

Outputs, per chromosome:
    - n_sequences, sequence_length, n_valid_sites
    - total_pairwise_differences   (sum_{i<j} k_ij -- raw count, exactly the
      numerator of the Nei & Li formula)
    - mean_pairwise_differences    (= total_pairwise_differences / (n(n-1)/2),
      i.e. pi exactly as written in the formula -- NOT divided by sequence
      length)
    - pi_per_site                 (mean_pairwise_differences / n_valid_sites,
      i.e. nucleotide diversity per site -- the conventional way pi is
      reported and compared across regions of different length)

A genome-wide summary (pooling all chromosomes, assuming the same n sequences
across all of them) is also reported.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
from Bio import SeqIO

VALID_BASES = [b"A", b"C", b"G", b"T"]


def load_chromosome_matrix(path):
    """Load a chromosome fasta into an (n_sequences, length) numpy byte matrix."""
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No records found in {path}")

    n = len(records)
    length = len(records[0].seq)

    seq_bytes = bytearray()
    for record in records:
        s = str(record.seq).upper()
        if len(s) != length:
            raise ValueError(
                f"{path}: record {record.id} has length {len(s)}, "
                f"expected {length} (based on the first record). "
                f"All sequences in a chromosome file must be the same length."
            )
        seq_bytes.extend(s.encode("ascii"))

    arr = np.frombuffer(bytes(seq_bytes), dtype="S1").reshape(n, length)
    return arr


def compute_pi(arr):
    """Compute Nei & Li (1979) pairwise-difference stats for one chromosome matrix.

    arr: numpy array of shape (n_sequences, length), dtype 'S1'.
    Returns a dict of summary stats.
    """
    n_seqs, length = arr.shape

    # per-site counts of each valid base -> each entry shape (length,)
    counts = {b: (arr == b).sum(axis=0) for b in VALID_BASES}

    valid_n = sum(counts.values())                                # called bases per site
    same_pairs = sum(c * (c - 1) // 2 for c in counts.values())   # same-allele pairs per site
    total_pairs = valid_n * (valid_n - 1) // 2                    # all possible pairs per site
    diff_pairs = total_pairs - same_pairs                         # differing pairs per site

    # only sites where at least 2 sequences have a called base contribute
    usable = valid_n >= 2
    n_valid_sites = int(usable.sum())

    total_pairwise_differences = int(diff_pairs[usable].sum())

    # pi exactly as written in the Nei & Li formula: sum(k_ij) / (n(n-1)/2),
    # using the full sample size n_seqs (assumes essentially-complete data --
    # sites with missing calls still contribute what differences they can,
    # they just don't inflate the n(n-1)/2 denominator here)
    n_pairs_full = n_seqs * (n_seqs - 1) // 2
    mean_pairwise_differences = (
        total_pairwise_differences / n_pairs_full if n_pairs_full else float("nan")
    )

    # nucleotide diversity per site (the conventional way pi is reported)
    pi_per_site = (
        mean_pairwise_differences / n_valid_sites if n_valid_sites else float("nan")
    )

    return {
        "n_sequences": n_seqs,
        "sequence_length": length,
        "n_valid_sites": n_valid_sites,
        "total_pairwise_differences": total_pairwise_differences,
        "mean_pairwise_differences": mean_pairwise_differences,
        "pi_per_site": pi_per_site,
    }


def chrom_sort_key(path):
    """Sort chr1, chr2, ..., chr10 numerically instead of lexicographically."""
    digits = re.findall(r"\d+", os.path.basename(path))
    return int(digits[0]) if digits else 0


def main():
    parser = argparse.ArgumentParser(
        description="Compute Nei & Li (1979) nucleotide diversity (pi) per chromosome."
    )
    parser.add_argument(
        "--chrom-dir",
        required=True,
        help="Directory containing chr1.fasta ... chrN.fasta",
    )
    parser.add_argument(
        "--pattern",
        default="chr*.fasta",
        help="Glob pattern (relative to --chrom-dir) for chromosome fasta files (default: chr*.fasta)",
    )
    parser.add_argument(
        "--out",
        default="diversity_summary.tsv",
        help="Output TSV path (default: diversity_summary.tsv)",
    )
    args = parser.parse_args()

    paths = sorted(
        glob.glob(os.path.join(args.chrom_dir, args.pattern)),
        key=chrom_sort_key,
    )
    if not paths:
        print(f"ERROR: no files matched {args.pattern} in {args.chrom_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(paths)} chromosome files: {[os.path.basename(p) for p in paths]}")

    rows = []
    grand_total_diff = 0
    grand_total_valid_sites = 0
    n_pairs_reference = None

    for path in paths:
        chrom_name = os.path.splitext(os.path.basename(path))[0]
        print(f"Processing {chrom_name} ...")

        arr = load_chromosome_matrix(path)
        stats = compute_pi(arr)
        stats["chromosome"] = chrom_name
        rows.append(stats)

        n_pairs_full = stats["n_sequences"] * (stats["n_sequences"] - 1) // 2
        if n_pairs_reference is None:
            n_pairs_reference = n_pairs_full
        elif n_pairs_full != n_pairs_reference:
            print(
                f"WARNING: {chrom_name} has {stats['n_sequences']} sequences, which "
                f"differs from an earlier chromosome. The genome-wide summary assumes "
                f"a consistent sample size across all chromosomes -- treat the "
                f"genome-wide row with caution if this warning appears.",
                file=sys.stderr,
            )

        grand_total_diff += stats["total_pairwise_differences"]
        grand_total_valid_sites += stats["n_valid_sites"]

        print(
            f"  {chrom_name}: n={stats['n_sequences']} L={stats['sequence_length']} "
            f"valid_sites={stats['n_valid_sites']} "
            f"mean_pairwise_diff={stats['mean_pairwise_differences']:.4f} "
            f"pi_per_site={stats['pi_per_site']:.6g}"
        )

    genome_mean_pairwise = (
        grand_total_diff / n_pairs_reference if n_pairs_reference else float("nan")
    )
    genome_pi_per_site = (
        genome_mean_pairwise / grand_total_valid_sites if grand_total_valid_sites else float("nan")
    )

    with open(args.out, "w") as f:
        f.write(
            "chromosome\tn_sequences\tsequence_length\tn_valid_sites\t"
            "total_pairwise_differences\tmean_pairwise_differences\tpi_per_site\n"
        )
        for stats in rows:
            f.write(
                f"{stats['chromosome']}\t{stats['n_sequences']}\t{stats['sequence_length']}\t"
                f"{stats['n_valid_sites']}\t{stats['total_pairwise_differences']}\t"
                f"{stats['mean_pairwise_differences']:.6f}\t{stats['pi_per_site']:.8f}\n"
            )
        f.write(
            f"GENOME_WIDE\tNA\tNA\t{grand_total_valid_sites}\t{grand_total_diff}\t"
            f"{genome_mean_pairwise:.6f}\t{genome_pi_per_site:.8f}\n"
        )

    print(f"\nDone. Wrote per-chromosome and genome-wide diversity to: {args.out}")
    print(f"Genome-wide pi (per site): {genome_pi_per_site:.8f}")


if __name__ == "__main__":
    main()
