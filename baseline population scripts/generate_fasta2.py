#!/usr/bin/env python3
"""
popgenART: Generate FASTA sequences per sample from a reference + SNP positions + SNP calls.

Fast replacement for the bash version:
  - Reference and SNP-position parsing happen ONCE (not per sample/chromosome).
  - Chromosome order is preserved from indices_SNP (fixes the offset-assignment
    bug caused by bash associative-array iteration being unordered).
  - Sequence mutation uses bytearray (O(1) per SNP) instead of string rebuilding.
  - Output is written incrementally, one sample at a time, so memory stays flat
    regardless of chunk size.

Usage:
  python3 generate_fasta.py -i indices_SNP.txt -r ref_seq.fa -s SNP_chunk.csv -o out_chunk.fasta
"""
import argparse
import sys


def parse_indices(path, one_based=False):
    """Parse indices_SNP into (chrom_order, {chrom: [positions...]}).
    chrom_order preserves the order chromosomes appear in the file -- this
    must match the order the SNP_csv Sequence column was built in.
    If one_based is True, every position is shifted down by 1 to convert
    to Python's 0-based indexing."""
    order = []
    positions = {}
    cur = None
    shift = 1 if one_based else 0
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("Chrom"):
                cur = line.split()[2]  # "Chrom ---- 1" -> "1"
                order.append(cur)
                positions[cur] = []
            else:
                positions[cur].append(int(line.strip()) - shift)
    return order, positions


def parse_ref(path):
    """Parse reference FASTA into {chrom_number_str: sequence_string}.
    Expects headers like '>Chromosome 1'; uses the last whitespace-delimited
    token as the chromosome key (same convention as the original awk match)."""
    seqs = {}
    cur = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(chunks)
                cur = line.split()[-1]
                chunks = []
            else:
                chunks.append(line)
        if cur is not None:
            seqs[cur] = "".join(chunks)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--indices", required=True, help="indices_SNP file")
    ap.add_argument("-r", "--ref", required=True, help="reference FASTA")
    ap.add_argument("-s", "--csv", required=True, help="SNP CSV chunk (with header row)")
    ap.add_argument("-o", "--out", required=True, help="output FASTA path for this chunk")
    ap.add_argument(
        "--one-based",
        action="store_true",
        help="treat positions in indices_SNP as 1-based genomic coordinates "
             "(shifts every position down by 1 for Python's 0-based indexing)",
    )
    args = ap.parse_args()

    print("Parsing SNP indices...", file=sys.stderr)
    chr_order, positions = parse_indices(args.indices, one_based=args.one_based)

    print(f"Parsing reference sequences ({len(chr_order)} chromosomes)...", file=sys.stderr)
    ref_seqs = parse_ref(args.ref)
    missing = [c for c in chr_order if c not in ref_seqs]
    if missing:
        sys.exit(f"ERROR: chromosomes in indices file but missing from reference: {missing[:10]}")
    ref_bytes = {c: ref_seqs[c].encode("ascii") for c in chr_order}

    print("Validating SNP positions against reference lengths...", file=sys.stderr)
    problems = []
    for chr_num in chr_order:
        ref_len = len(ref_bytes[chr_num])
        for pos in positions[chr_num]:
            if pos < 0 or pos >= ref_len:
                problems.append((chr_num, pos, ref_len))
    if problems:
        print(f"ERROR: {len(problems)} SNP position(s) fall outside their "
              f"reference sequence bounds. First 10 shown:", file=sys.stderr)
        for chr_num, pos, ref_len in problems[:10]:
            print(f"  chr{chr_num}: position {pos} (valid range 0-{ref_len - 1}, "
                  f"ref length {ref_len})", file=sys.stderr)
        print(
            "If your indices_SNP positions are 1-based genomic coordinates, "
            "rerun with --one-based. Otherwise check the reference lengths "
            "and indices file for a mismatch (e.g. wrong reference version).",
            file=sys.stderr,
        )
        sys.exit(1)

    print("Processing samples...", file=sys.stderr)
    n_samples = 0
    with open(args.csv) as f, open(args.out, "w") as out:
        header = f.readline()  # discard header
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            SampleID, _col2, Sequence = line.split(",", 2)
            seq_bytes = Sequence.encode("ascii")
            offset = 0
            for chr_num in chr_order:
                arr = bytearray(ref_bytes[chr_num])  # fresh mutable copy
                for pos in positions[chr_num]:
                    arr[pos] = seq_bytes[offset]
                    offset += 1
                out.write(f">{SampleID}chr{chr_num}\n")
                out.write(arr.decode("ascii"))
                out.write("\n")
            n_samples += 1
            if offset != len(seq_bytes):
                print(
                    f"WARNING: sample {SampleID} consumed {offset} of "
                    f"{len(seq_bytes)} SNP characters (length mismatch).",
                    file=sys.stderr,
                )

    print(f"Done: {n_samples} samples written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()