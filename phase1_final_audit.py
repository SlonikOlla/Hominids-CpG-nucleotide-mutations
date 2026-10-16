#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path("results/phase1_genome")
AGG = ROOT / "aggregated"
CHROMS = [f"chr{i}" for i in range(1, 23)]
LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
THRESHOLDS = ["0.90", "0.95", "0.99"]
CONTEXTS = ["CG", "CHG", "CHH"]
COUNTS = ("opportunities", "C>T", "C>A",
          "C>G", "total_substitutions")
EXPECTED_MOTIFS = {"CG": 64, "CHG": 48, "CHH": 144}

def read(path):
    if not path.is_file():
        raise AssertionError(f"Missing: {path}")
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def key(row):
    return (row["lineage"], row["posterior_threshold"],
            row["context"], row["ancestral_5mer"])

def validate(row):
    assert row["lineage"] in LINEAGES
    assert row["posterior_threshold"] in THRESHOLDS
    assert row["context"] in CONTEXTS
    m = row["ancestral_5mer"]
    assert len(m) == 5 and m[2] == "C"
    assert set(m) <= set("ACGT")
    context = ("CG" if m[3] == "G" else
               "CHG" if m[4] == "G" else "CHH")
    assert context == row["context"]
    c = [int(row[x]) for x in COUNTS]
    assert all(x >= 0 for x in c)
    assert c[4] == sum(c[1:4])
    assert c[4] <= c[0]
    n = c[0]
    if n:
        assert math.isclose(
            float(row["rate_per_million"]),
            c[4] * 1e6 / n, abs_tol=0.000002)
        assert math.isclose(
            float(row["CT_rate_per_million"]),
            c[1] * 1e6 / n, abs_tol=0.000002)

def sums():
    return {x: 0 for x in COUNTS}

chrom_sum = defaultdict(sums)
context_sum = defaultdict(sums)

for chrom in CHROMS:
    path = ROOT / f"{chrom}_formal_5mer_STRICT_summary.production.tsv"
    qc = ROOT / f"{chrom}_formal_5mer_STRICT_QC.production.tsv"
    assert qc.is_file(), f"Missing QC: {qc}"
    rows = read(path)
    assert len(rows) == 3072, (chrom, len(rows))
    seen = set()
    motif_sets = defaultdict(set)

    for row in rows:
        validate(row)
        k = key(row)
        assert k not in seen, (chrom, k)
        seen.add(k)
        lineage, threshold, context, motif = k
        motif_sets[(lineage, threshold, context)].add(motif)
        for field in COUNTS:
            value = int(row[field])
            chrom_sum[k][field] += value
            context_sum[k[:3]][field] += value

    for lineage in LINEAGES:
        for threshold in THRESHOLDS:
            for context in CONTEXTS:
                assert len(motif_sets[(lineage, threshold, context)]) == \
                    EXPECTED_MOTIFS[context]
    print("PASS:", chrom, "3072 rows")

motif_rows = read(AGG / "genome_exact_5mer_summary.tsv")
assert len(motif_rows) == 3072
seen = set()

for row in motif_rows:
    validate(row)
    k = key(row)
    assert k not in seen
    seen.add(k)
    assert all(int(row[f]) == chrom_sum[k][f] for f in COUNTS)

assert seen == set(chrom_sum)

context_rows = read(AGG / "genome_context_summary.tsv")
assert len(context_rows) == 36
seen_context = set()
observed = {}

for row in context_rows:
    k = (row["lineage"], row["posterior_threshold"], row["context"])
    assert k not in seen_context
    seen_context.add(k)
    assert all(int(row[f]) == context_sum[k][f] for f in COUNTS)
    n = int(row["opportunities"])
    ct = int(row["C>T"])
    total = int(row["total_substitutions"])
    assert n > 0
    assert math.isclose(float(row["CT_rate_per_million"]),
                        ct * 1e6 / n, abs_tol=0.000002)
    assert math.isclose(float(row["rate_per_million"]),
                        total * 1e6 / n, abs_tol=0.000002)
    observed[k] = (n, ct)

assert len(seen_context) == 36

for lineage in LINEAGES:
    for context in CONTEXTS:
        ns = [observed[(lineage, t, context)][0]
              for t in THRESHOLDS]
        assert ns[0] >= ns[1] >= ns[2], (lineage, context, ns)

ratio_rows = read(AGG / "genome_context_ratios.tsv")
assert len(ratio_rows) == 24
seen_ratios = set()

for row in ratio_rows:
    lineage = row["lineage"]
    threshold = row["posterior_threshold"]
    comparison = row["comparison"]
    k = (lineage, threshold, comparison)
    assert k not in seen_ratios
    seen_ratios.add(k)
    assert comparison in ("CG/CHG", "CG/CHH")
    denominator = comparison.split("/")[1]

    def get_rate(context, mutation):
        c = context_sum[(lineage, threshold, context)]
        return c[mutation] / c["opportunities"]

    expected_ct = get_rate("CG", "C>T") / \
                  get_rate(denominator, "C>T")
    expected_total = get_rate("CG", "total_substitutions") / \
                     get_rate(denominator, "total_substitutions")
    assert math.isclose(float(row["CT_rate_ratio"]),
                        expected_ct, abs_tol=0.000002)
    assert math.isclose(float(row["total_rate_ratio"]),
                        expected_total, abs_tol=0.000002)

print("PASS: all aggregated 5-mer counts")
print("PASS: all aggregated context counts")
print("PASS: all rate calculations")
print("PASS: all context ratios")
print("PASS: posterior opportunity monotonicity")
print("PASS: all 22 chromosome summary and QC files present")
print("FINAL AUDIT: PASS")
