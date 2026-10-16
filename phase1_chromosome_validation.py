#!/usr/bin/env python3

import csv
import random
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path("results/phase1_genome/aggregated")
SOURCE = ROOT / "chromosome_context_summary.tsv"

LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
CHROMS = [f"chr{i}" for i in range(1, 23)]
CONTEXTS = ["CG", "CHG", "CHH"]
B = 10000

data = {}

with open(SOURCE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["posterior_threshold"] != "0.95":
            continue

        key = (
            row["chromosome"],
            row["lineage"],
            row["context"]
        )

        if key in data:
            raise ValueError(f"Duplicate row: {key}")

        data[key] = {
            "opportunities": int(row["opportunities"]),
            "CT": int(row["C>T"])
        }

for chrom in CHROMS:
    for lineage in LINEAGES:
        for context in CONTEXTS:
            assert (chrom, lineage, context) in data

def rate(counts):
    if counts["opportunities"] == 0:
        raise ValueError("Zero opportunities")
    return counts["CT"] / counts["opportunities"]

def ratio(chrom, lineage, denominator):
    return (
        rate(data[(chrom, lineage, "CG")])
        / rate(data[(chrom, lineage, denominator)])
    )

def pooled_ratio(chroms, lineage, denominator):
    def pooled_rate(context):
        ct = sum(
            data[(c, lineage, context)]["CT"]
            for c in chroms
        )
        opportunities = sum(
            data[(c, lineage, context)]["opportunities"]
            for c in chroms
        )
        return ct / opportunities

    return pooled_rate("CG") / pooled_rate(denominator)

def quantile(values, q):
    values = sorted(values)
    x = (len(values) - 1) * q
    lo = int(x)
    hi = min(lo + 1, len(values) - 1)
    return values[lo] + (x - lo) * (values[hi] - values[lo])

rng = random.Random(20261008)

chrom_rows = []
summary_rows = []

for lineage in LINEAGES:
    for denominator in ("CHG", "CHH"):
        ratios = []

        for chrom in CHROMS:
            r = ratio(chrom, lineage, denominator)
            ratios.append(r)

            chrom_rows.append({
                "chromosome": chrom,
                "lineage": lineage,
                "comparison": f"CG/{denominator}",
                "CT_rate_ratio": f"{r:.6f}"
            })

        bootstrap = []

        for _ in range(B):
            sampled = rng.choices(CHROMS, k=len(CHROMS))
            bootstrap.append(
                pooled_ratio(sampled, lineage, denominator)
            )

        observed = pooled_ratio(CHROMS, lineage, denominator)

        summary_rows.append({
            "lineage": lineage,
            "comparison": f"CG/{denominator}",
            "pooled_ratio": f"{observed:.6f}",
            "chromosome_median": f"{statistics.median(ratios):.6f}",
            "chromosome_min": f"{min(ratios):.6f}",
            "chromosome_max": f"{max(ratios):.6f}",
            "bootstrap_CI_low": f"{quantile(bootstrap, 0.025):.6f}",
            "bootstrap_CI_high": f"{quantile(bootstrap, 0.975):.6f}",
            "chromosomes_with_CG_enrichment": sum(
                r > 1 for r in ratios
            )
        })

def write(name, rows):
    path = ROOT / name
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    print("WROTE:", path)

write("chromosome_CT_ratios.tsv", chrom_rows)
write("chromosome_validation_summary.tsv", summary_rows)

print("\n===== CHROMOSOME VALIDATION =====")
for row in summary_rows:
    print(
        row["lineage"],
        row["comparison"],
        "pooled:", row["pooled_ratio"],
        "95% bootstrap interval:",
        row["bootstrap_CI_low"],
        row["bootstrap_CI_high"],
        "enriched chromosomes:",
        row["chromosomes_with_CG_enrichment"], "/22"
    )

print("\nSTATUS: PASS")
