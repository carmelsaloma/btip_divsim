import random
import sys

# ---------------------------
# Config (seed can be passed as an argument from the bash script)
# ---------------------------
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
random.seed(seed)

N_CHROMOSOMES  = 10
LOCI_PER_CHR   = 150
MIN_NONNEUTRAL = 8

NEUTRAL_RANGE    = range(1, 1401)     # 1 - 1400
NONNEUTRAL_RANGE = range(1401, 1501)  # 1401 - 1500

OUT_TXT = "/home/projects/btip26_divsim/mapping150/chromosome_loci_reference.txt"

# ---------------------------
# Build and shuffle pools
# ---------------------------
neutral    = list(NEUTRAL_RANGE)
nonneutral = list(NONNEUTRAL_RANGE)
random.shuffle(neutral)
random.shuffle(nonneutral)

chromosomes = {}

# step 1: guarantee minimum 8 non-neutral loci per chromosome
for c in range(1, N_CHROMOSOMES + 1):
    chromosomes[f"chr{c}"] = nonneutral[:MIN_NONNEUTRAL]
    nonneutral = nonneutral[MIN_NONNEUTRAL:]

# step 2: mix leftover non-neutral with all neutral loci, fill each chromosome to 150
leftover = nonneutral + neutral
random.shuffle(leftover)

for c in range(1, N_CHROMOSOMES + 1):
    key = f"chr{c}"
    slots_left = LOCI_PER_CHR - len(chromosomes[key])
    chromosomes[key] += leftover[:slots_left]
    leftover = leftover[slots_left:]
    random.shuffle(chromosomes[key])

# ---------------------------
# Write reference .txt file
# ---------------------------
with open(OUT_TXT, "w") as f:
    f.write(f"Chromosome loci assignment reference (seed={seed})\n")
    f.write("=" * 55 + "\n\n")
    for c in range(1, N_CHROMOSOMES + 1):
        key = f"chr{c}"
        loci = chromosomes[key]
        n_nonneutral = sum(1 for x in loci if x > 1400)
        n_neutral = sum(1 for x in loci if x <= 1400)
        f.write(f"{key}: total={len(loci)}  neutral={n_neutral}  nonneutral={n_nonneutral}\n")
        f.write(", ".join(str(x) for x in sorted(loci)) + "\n\n")

    total_used = sum(len(v) for v in chromosomes.values())
    f.write(f"Total loci used: {total_used} / 1500\n")

print(f"Reference file written to: {OUT_TXT}")

# also print a quick summary to console (visible in bash output)
for c in range(1, N_CHROMOSOMES + 1):
    key = f"chr{c}"
    loci = chromosomes[key]
    n_nonneutral = sum(1 for x in loci if x > 1400)
    n_neutral = sum(1 for x in loci if x <= 1400)
    print(f"{key}: total={len(loci)}  neutral={n_neutral}  nonneutral={n_nonneutral}")
