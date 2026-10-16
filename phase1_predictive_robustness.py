#!/usr/bin/env python3

import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(
    "results/phase1_strengthening/predictive_comparison"
)

INPUT = ROOT / "heldout_chromosome_predictions.tsv"

LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
CONTEXTS = ["CG", "CHG", "CHH"]
CHROMS = [f"chr{i}" for i in range(1, 23)]

rows = []

with open(INPUT, newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append({
            "chromosome": r["chromosome"],
            "lineage": r["lineage"],
            "context": r["context"],
            "motif": r["motif"],
            "n": int(r["opportunities"]),
            "k": int(r["CT"]),
            "ll3": float(r["LL_3mer"]),
            "ll5": float(r["LL_5mer"])
        })

assert len(rows) == 22 * 4 * 256, (
    f"Unexpected number of rows: {len(rows)}"
)

def saturated_ll(k, n):
    """Binomial saturated log likelihood, omitting common constants."""
    if k == 0 or k == n:
        return 0.0

    p = k / n
    return k * math.log(p) + (n-k) * math.log1p(-p)

def summarize(subset):
    n = sum(r["n"] for r in subset)

    ll3 = sum(r["ll3"] for r in subset)
    ll5 = sum(r["ll5"] for r in subset)

    llsat = sum(
        saturated_ll(r["k"], r["n"])
        for r in subset
    )

    gain = ll5 - ll3

    dev3 = 2 * (llsat - ll3)
    dev5 = 2 * (llsat - ll5)

    by_chrom = defaultdict(float)

    for r in subset:
        by_chrom[r["chromosome"]] += (
            r["ll5"] - r["ll3"]
        )

    return {
        "opportunities": n,
        "LL_gain": gain,
        "LL_gain_per_million": gain / n * 1e6,
        "deviance_3mer": dev3,
        "deviance_5mer": dev5,
        "deviance_reduction": dev3 - dev5,
        "percent_deviance_reduction": (
            100 * (dev3 - dev5) / dev3
            if dev3 > 0 else float("nan")
        ),
        "positive_chromosomes": sum(
            by_chrom[c] > 0 for c in CHROMS
        )
    }

output = []

for lineage in LINEAGES:
    for context in CONTEXTS:
        for exclude in (False, True):
            subset = [
                r for r in rows
                if r["lineage"] == lineage
                and r["context"] == context
                and (not exclude or r["motif"] != "TTCTT")
            ]

            result = summarize(subset)

            output.append({
                "lineage": lineage,
                "context": context,
                "exclude_TTCTT": exclude,
                **result
            })

    for exclude in (False, True):
        subset = [
            r for r in rows
            if r["lineage"] == lineage
            and (not exclude or r["motif"] != "TTCTT")
        ]

        output.append({
            "lineage": lineage,
            "context": "ALL",
            "exclude_TTCTT": exclude,
            **summarize(subset)
        })

path = ROOT / "predictive_robustness_summary.tsv"

with open(path, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(output[0]),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(output)

print("PREDICTIVE ROBUSTNESS ANALYSIS")
print("Input rows:", len(rows))

print("\nPER-MILLION IMPROVEMENT, TTCTT INCLUDED")
print(
    "lineage   context   LL_gain_per_million"
    "   deviance_reduction_%   positive_chromosomes"
)

for r in output:
    if r["exclude_TTCTT"]:
        continue

    print(
        f"{r['lineage']:8s} "
        f"{r['context']:7s} "
        f"{r['LL_gain_per_million']:19.3f} "
        f"{r['percent_deviance_reduction']:22.3f} "
        f"{r['positive_chromosomes']:>3}/22"
    )

print("\nTTCTT EXCLUSION TEST")
print(
    "lineage   context  "
    "full_gain    excluded_gain   retained_%"
    "   positive_chromosomes"
)

for lineage in LINEAGES:
    for context in ("CHH", "ALL"):
        full = next(
            r for r in output
            if r["lineage"] == lineage
            and r["context"] == context
            and not r["exclude_TTCTT"]
        )

        excluded = next(
            r for r in output
            if r["lineage"] == lineage
            and r["context"] == context
            and r["exclude_TTCTT"]
        )

        retained = (
            100 * excluded["LL_gain"] / full["LL_gain"]
        )

        print(
            f"{lineage:8s} "
            f"{context:7s} "
            f"{full['LL_gain']:10.1f} "
            f"{excluded['LL_gain']:14.1f} "
            f"{retained:10.2f} "
            f"{excluded['positive_chromosomes']:>3}/22"
        )

print("\nSaved:", path)
