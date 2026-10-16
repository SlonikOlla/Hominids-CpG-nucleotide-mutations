#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path("results/phase1_genome")
OUT = Path("results/phase1_strengthening/predictive_comparison")
OUT.mkdir(parents=True, exist_ok=True)

LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
CONTEXTS = ["CG", "CHG", "CHH"]
CHROMS = [f"chr{i}" for i in range(1, 23)]
PP = 0.95

data = defaultdict(dict)

for chrom in CHROMS:
    path = ROOT / f"{chrom}_formal_5mer_STRICT_summary.production.tsv"
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if abs(float(r["posterior_threshold"]) - PP) > 1e-8:
                continue
            key = (r["lineage"], r["ancestral_5mer"])
            assert key not in data[chrom]
            data[chrom][key] = (
                r["context"],
                int(r["opportunities"]),
                int(r["C>T"])
            )

motifs = sorted({m for _, m in data[CHROMS[0]]})
assert len(motifs) == 256

for chrom in CHROMS:
    assert len(data[chrom]) == 4 * 256

def probability(successes, opportunities):
    # Jeffreys smoothing avoids zero-probability predictions.
    return (successes + 0.5) / (opportunities + 1.0)

def log_likelihood(successes, opportunities, p):
    return (successes * math.log(p) +
            (opportunities - successes) * math.log1p(-p))

rows = []

for heldout in CHROMS:
    train3 = defaultdict(lambda: [0, 0])
    train5 = defaultdict(lambda: [0, 0])

    for chrom in CHROMS:
        if chrom == heldout:
            continue
        for (lineage, motif), (context, n, k) in data[chrom].items():
            key3 = (lineage, motif[1:4])
            key5 = (lineage, motif)
            for store, key in ((train3, key3), (train5, key5)):
                store[key][0] += n
                store[key][1] += k

    for (lineage, motif), (context, n, k) in data[heldout].items():
        n3, k3 = train3[(lineage, motif[1:4])]
        n5, k5 = train5[(lineage, motif)]
        p3 = probability(k3, n3)
        p5 = probability(k5, n5)

        ll3 = log_likelihood(k, n, p3)
        ll5 = log_likelihood(k, n, p5)

        rows.append({
            "chromosome": heldout,
            "lineage": lineage,
            "context": context,
            "motif": motif,
            "opportunities": n,
            "CT": k,
            "LL_3mer": ll3,
            "LL_5mer": ll5,
            "LL_gain_5mer": ll5 - ll3
        })

path = OUT / "heldout_chromosome_predictions.tsv"

with open(path, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("LEAVE-ONE-CHROMOSOME-OUT COMPARISON")
print("Posterior threshold:", PP)
print("Held-out chromosomes:", len(CHROMS))
print("3-mer categories:", len({m[1:4] for m in motifs}))
print("5-mer categories:", len(motifs))
print()
print("LINEAGE  CONTEXT  LL_GAIN_5MER  DEVIANCE_REDUCTION  POSITIVE_CHROMS")

for lineage in LINEAGES:
    for context in CONTEXTS:
        subset = [
            r for r in rows
            if r["lineage"] == lineage and r["context"] == context
        ]
        gain = sum(r["LL_gain_5mer"] for r in subset)
        positive = sum(
            sum(r["LL_gain_5mer"] for r in subset
                if r["chromosome"] == chrom) > 0
            for chrom in CHROMS
        )
        print(
            f"{lineage:8s} {context:7s} "
            f"{gain:14.2f} {2*gain:19.2f} "
            f"{positive:>3}/22"
        )

print()
print("OVERALL BY LINEAGE")
for lineage in LINEAGES:
    subset = [r for r in rows if r["lineage"] == lineage]
    gain = sum(r["LL_gain_5mer"] for r in subset)
    positive = sum(
        sum(r["LL_gain_5mer"] for r in subset
            if r["chromosome"] == chrom) > 0
        for chrom in CHROMS
    )
    print(
        f"{lineage:8s} LL_gain={gain:.2f} "
        f"deviance_reduction={2*gain:.2f} "
        f"positive_chromosomes={positive}/22"
    )

print("\nOutput:", path)
