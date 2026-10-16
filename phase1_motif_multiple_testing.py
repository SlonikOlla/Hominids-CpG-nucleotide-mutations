#!/usr/bin/env python3

import csv
import itertools
from pathlib import Path
import numpy as np

ROOT = Path("results/phase1_genome")
OUT = Path("results/phase1_strengthening/multiple_testing")
OUT.mkdir(parents=True, exist_ok=True)

LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
CONTEXTS = ["CG", "CHG", "CHH"]
CHROMS = [f"chr{i}" for i in range(1, 23)]

B = 10000
SEED = 20261008
THRESHOLD = 0.95

# Read all 22 chromosome summaries.
records = {}
motif_context = {}

for chrom in CHROMS:
    path = ROOT / f"{chrom}_formal_5mer_STRICT_summary.production.tsv"

    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if abs(float(r["posterior_threshold"]) - THRESHOLD) > 1e-8:
                continue

            motif = r["ancestral_5mer"]
            context = r["context"]
            lineage = r["lineage"]

            key = (chrom, lineage, motif)

            if key in records:
                raise ValueError(f"Duplicate record: {key}")

            records[key] = (
                int(r["opportunities"]),
                int(r["C>T"])
            )
            motif_context[motif] = context

motifs = sorted(motif_context)
assert len(motifs) == 256

opp = np.zeros((22, 4, 256))
ct = np.zeros_like(opp)

for c, chrom in enumerate(CHROMS):
    for l, lineage in enumerate(LINEAGES):
        for m, motif in enumerate(motifs):
            key = (chrom, lineage, motif)
            if key not in records:
                raise ValueError(f"Missing: {key}")

            opp[c, l, m], ct[c, l, m] = records[key]

masks = {
    context: np.array([
        motif_context[m] == context for m in motifs
    ])
    for context in CONTEXTS
}

def effects(o, s):
    """Context-normalized log2 C>T motif effects."""
    result = np.full((4, 256), np.nan)

    for context, mask in masks.items():
        base_o = o[:, mask].sum(axis=1)
        base_s = s[:, mask].sum(axis=1)

        base = base_s / base_o

        with np.errstate(divide="ignore", invalid="ignore"):
            rates = s[:, mask] / o[:, mask]
            result[:, mask] = np.log2(rates / base[:, None])

    return result

pairs = list(itertools.combinations(range(4), 2))

def contrasts(o, s):
    e = effects(o, s)
    return np.stack([
        e[a] - e[b] for a, b in pairs
    ])

observed = contrasts(opp.sum(axis=0), ct.sum(axis=0))

assert np.all(np.isfinite(observed)), \
    "Nonfinite observed contrasts"

rng = np.random.default_rng(SEED)

boot = np.empty((B, len(pairs), len(motifs)))

for b in range(B):
    sampled = rng.integers(0, 22, size=22)
    o = opp[sampled].sum(axis=0)
    s = ct[sampled].sum(axis=0)
    boot[b] = contrasts(o, s)

assert np.all(np.isfinite(boot)), \
    "Nonfinite bootstrap contrasts"

# Pointwise uncertainty.
point_lo = np.quantile(boot, 0.025, axis=0)
point_hi = np.quantile(boot, 0.975, axis=0)

# Studentized maximum-error bootstrap.
# Center each bootstrap contrast on its observed value.
errors = boot - observed[None, :, :]

se = boot.std(axis=0, ddof=1)
assert np.all(se > 0), "Zero bootstrap standard error"

standardized = np.abs(errors / se[None, :, :])

# One maximum across all 6 x 256 comparisons per replicate.
max_t = standardized.reshape(B, -1).max(axis=1)
critical = np.quantile(max_t, 0.95)

sim_lo = observed - critical * se
sim_hi = observed + critical * se

# Multiplicity-adjusted max-T bootstrap p values.
abs_t = np.abs(observed / se)
adjusted_p = np.mean(
    max_t[:, None, None] >= abs_t[None, :, :],
    axis=0
)

rows = []

for p, (a, b) in enumerate(pairs):
    for m, motif in enumerate(motifs):
        rows.append({
            "context": motif_context[motif],
            "motif": motif,
            "lineage_1": LINEAGES[a],
            "lineage_2": LINEAGES[b],
            "delta_log2": float(observed[p, m]),
            "bootstrap_SE": float(se[p, m]),
            "pointwise_CI_low": float(point_lo[p, m]),
            "pointwise_CI_high": float(point_hi[p, m]),
            "simultaneous_CI_low": float(sim_lo[p, m]),
            "simultaneous_CI_high": float(sim_hi[p, m]),
            "maxT_adjusted_p": float(adjusted_p[p, m]),
            "simultaneous_CI_excludes_zero": bool(
                sim_lo[p, m] > 0 or sim_hi[p, m] < 0
            )
        })

rows.sort(key=lambda r: r["maxT_adjusted_p"])

path = OUT / "motif_pairwise_simultaneous_inference.tsv"

with open(path, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

significant = [
    r for r in rows
    if r["simultaneous_CI_excludes_zero"]
]

print("MULTIPLE-TESTING ANALYSIS COMPLETE")
print("Chromosomes:", len(CHROMS))
print("Motifs:", len(motifs))
print("Lineage pairs:", len(pairs))
print("Total comparisons:", len(rows))
print("Bootstrap replicates:", B)
print("Simultaneous critical value:", round(critical, 4))
print("Simultaneous intervals excluding zero:", len(significant))

print("\nTOP 20 BY ABSOLUTE EFFECT")
for r in sorted(
    rows,
    key=lambda x: abs(x["delta_log2"]),
    reverse=True
)[:20]:
    print(
        r["context"], r["motif"],
        r["lineage_1"], r["lineage_2"],
        f"delta={r['delta_log2']:.4f}",
        f"simCI=[{r['simultaneous_CI_low']:.4f},"
        f"{r['simultaneous_CI_high']:.4f}]",
        f"adjusted_p={r['maxT_adjusted_p']:.4f}"
    )

print("\nTTCTT COMPARISONS")
for r in rows:
    if r["motif"] == "TTCTT":
        print(
            r["lineage_1"], "vs", r["lineage_2"],
            f"delta={r['delta_log2']:.4f}",
            f"simCI=[{r['simultaneous_CI_low']:.4f},"
            f"{r['simultaneous_CI_high']:.4f}]",
            f"adjusted_p={r['maxT_adjusted_p']:.4f}"
        )

print("\nOutput:", path)
