#!/usr/bin/env python3

import csv
import itertools
import math
from pathlib import Path

ROOT = Path("results/phase1_genome/aggregated")
SOURCE = ROOT / "genome_5mer_validation.tsv"
OUT = ROOT / "motif_lineage_correlations.tsv"

LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
CONTEXTS = ["CG", "CHG", "CHH"]

data = {}

with open(SOURCE) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        key = (
            row["context"],
            row["lineage"],
            row["ancestral_5mer"]
        )
        if key in data:
            raise ValueError(f"Duplicate motif: {key}")
        data[key] = float(row["CT_rate_per_million"])

def ranks(values):
    ordered = sorted(enumerate(values), key=lambda x: x[1])
    result = [0.0] * len(values)
    i = 0

    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][1] == ordered[i][1]:
            j += 1
        average_rank = (i + 1 + j) / 2
        for k in range(i, j):
            result[ordered[k][0]] = average_rank
        i = j

    return result

def pearson(x, y):
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    numerator = sum((a-mx)*(b-my) for a,b in zip(x,y))
    denominator = math.sqrt(
        sum((a-mx)**2 for a in x) *
        sum((b-my)**2 for b in y)
    )
    return numerator / denominator

rows = []

for context in CONTEXTS:
    motifs = sorted({
        motif for ctx, lineage, motif in data
        if ctx == context
    })

    for a, b in itertools.combinations(LINEAGES, 2):
        x = [data[(context, a, motif)] for motif in motifs]
        y = [data[(context, b, motif)] for motif in motifs]

        rho = pearson(ranks(x), ranks(y))

        rows.append({
            "context": context,
            "lineage_1": a,
            "lineage_2": b,
            "motifs": len(motifs),
            "spearman_rho": f"{rho:.6f}"
        })

with open(OUT, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("===== CROSS-LINEAGE 5-MER CORRELATIONS =====")
for row in rows:
    print(
        row["context"],
        row["lineage_1"],
        row["lineage_2"],
        "rho =", row["spearman_rho"]
    )

print("\nWROTE:", OUT)
print("STATUS: PASS")
