#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

ROOT = Path("/project/def-vasilisa/hominid_methyl_evolution")
SOURCE = ROOT / "results/phase1_genome"
OUT = SOURCE / "aggregated"
OUT.mkdir(parents=True, exist_ok=True)

CHROMS = [f"chr{i}" for i in range(1, 23)]
LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
THRESHOLDS = ["0.90", "0.95", "0.99"]
CONTEXTS = ["CG", "CHG", "CHH"]
COUNT_FIELDS = [
    "opportunities", "C>T", "C>A", "C>G",
    "total_substitutions"
]

def accumulator():
    return {field: 0 for field in COUNT_FIELDS}

def add(target, row):
    for field in COUNT_FIELDS:
        target[field] += int(row[field])

def rates(row):
    n = row["opportunities"]
    if n == 0:
        return 0.0, 0.0
    return (
        row["total_substitutions"] * 1_000_000 / n,
        row["C>T"] * 1_000_000 / n
    )

def write_table(path, columns, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=columns, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    print("WROTE:", path)

genome_context = defaultdict(accumulator)
genome_5mer = defaultdict(accumulator)
chrom_context = defaultdict(accumulator)

print("===== READING 22 AUTOSOMES =====")

for chrom in CHROMS:
    path = SOURCE / (
        f"{chrom}_formal_5mer_STRICT_summary.production.tsv"
    )

    if not path.is_file():
        raise SystemExit(f"ERROR: missing {path}")

    nrows = 0

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            lineage = row["lineage"]
            threshold = row["posterior_threshold"]
            context = row["context"]
            motif = row["ancestral_5mer"]

            if lineage not in LINEAGES:
                raise ValueError((chrom, lineage))
            if threshold not in THRESHOLDS:
                raise ValueError((chrom, threshold))
            if context not in CONTEXTS:
                raise ValueError((chrom, context))

            if int(row["total_substitutions"]) != sum(
                int(row[k]) for k in ("C>T", "C>A", "C>G")
            ):
                raise ValueError(
                    f"Mutation count mismatch: {chrom} {row}"
                )

            add(
                genome_context[(lineage, threshold, context)],
                row
            )
            add(
                genome_5mer[
                    (lineage, threshold, context, motif)
                ],
                row
            )
            add(
                chrom_context[
                    (chrom, lineage, threshold, context)
                ],
                row
            )

            nrows += 1

    print(chrom, "rows:", nrows)

print("\n===== GENOME-WIDE CONTEXT SUMMARY =====")

columns = [
    "lineage", "posterior_threshold", "context",
    *COUNT_FIELDS, "rate_per_million",
    "CT_rate_per_million"
]

context_rows = []

for lineage in LINEAGES:
    for threshold in THRESHOLDS:
        for context in CONTEXTS:
            counts = genome_context[
                (lineage, threshold, context)
            ]
            rate, ct_rate = rates(counts)

            context_rows.append({
                "lineage": lineage,
                "posterior_threshold": threshold,
                "context": context,
                **counts,
                "rate_per_million": f"{rate:.6f}",
                "CT_rate_per_million": f"{ct_rate:.6f}"
            })

write_table(
    OUT / "genome_context_summary.tsv",
    columns, context_rows
)

print("\n===== EXACT 5-MER SUMMARY =====")

motif_columns = [
    "lineage", "posterior_threshold",
    "context", "ancestral_5mer",
    *COUNT_FIELDS,
    "rate_per_million", "CT_rate_per_million"
]

motif_rows = []

for key in sorted(genome_5mer):
    lineage, threshold, context, motif = key
    counts = genome_5mer[key]
    rate, ct_rate = rates(counts)

    motif_rows.append({
        "lineage": lineage,
        "posterior_threshold": threshold,
        "context": context,
        "ancestral_5mer": motif,
        **counts,
        "rate_per_million": f"{rate:.6f}",
        "CT_rate_per_million": f"{ct_rate:.6f}"
    })

write_table(
    OUT / "genome_exact_5mer_summary.tsv",
    motif_columns, motif_rows
)

print("\n===== CHROMOSOME SUMMARY =====")

chrom_columns = [
    "chromosome", *columns
]

chrom_rows = []

for chrom in CHROMS:
    for lineage in LINEAGES:
        for threshold in THRESHOLDS:
            for context in CONTEXTS:
                counts = chrom_context[
                    (chrom, lineage, threshold, context)
                ]
                rate, ct_rate = rates(counts)

                chrom_rows.append({
                    "chromosome": chrom,
                    "lineage": lineage,
                    "posterior_threshold": threshold,
                    "context": context,
                    **counts,
                    "rate_per_million": f"{rate:.6f}",
                    "CT_rate_per_million": f"{ct_rate:.6f}"
                })

write_table(
    OUT / "chromosome_context_summary.tsv",
    chrom_columns, chrom_rows
)

print("\n===== CONTEXT RATIOS =====")

ratio_columns = [
    "lineage", "posterior_threshold",
    "comparison", "total_rate_ratio", "CT_rate_ratio"
]

ratio_rows = []

for lineage in LINEAGES:
    for threshold in THRESHOLDS:
        r = {
            context: rates(
                genome_context[(lineage, threshold, context)]
            )
            for context in CONTEXTS
        }

        for denominator in ("CHG", "CHH"):
            total_ratio = (
                r["CG"][0] / r[denominator][0]
                if r[denominator][0] else float("nan")
            )
            ct_ratio = (
                r["CG"][1] / r[denominator][1]
                if r[denominator][1] else float("nan")
            )

            ratio_rows.append({
                "lineage": lineage,
                "posterior_threshold": threshold,
                "comparison": f"CG/{denominator}",
                "total_rate_ratio": f"{total_ratio:.6f}",
                "CT_rate_ratio": f"{ct_ratio:.6f}"
            })

write_table(
    OUT / "genome_context_ratios.tsv",
    ratio_columns, ratio_rows
)

print("\n===== PRIMARY RESULTS: PP >= 0.95 =====")

for lineage in LINEAGES:
    print("\n", lineage.upper())

    for context in CONTEXTS:
        c = genome_context[(lineage, "0.95", context)]
        rate, ct_rate = rates(c)

        print(
            context,
            "opportunities=", c["opportunities"],
            "C>T=", c["C>T"],
            "total=", c["total_substitutions"],
            "rate/M=", round(rate, 3),
            "CT/M=", round(ct_rate, 3)
        )

    for denominator in ("CHG", "CHH"):
        numerator = rates(
            genome_context[(lineage, "0.95", "CG")]
        )
        denominator_rates = rates(
            genome_context[(lineage, "0.95", denominator)]
        )

        print(
            f"CG/{denominator}:",
            "total ratio=",
            round(numerator[0] / denominator_rates[0], 3),
            "CT ratio=",
            round(numerator[1] / denominator_rates[1], 3)
        )

print("\nSTATUS: PASS")
