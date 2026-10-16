#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path("results/phase1_genome/aggregated")
SOURCE = ROOT / "genome_exact_5mer_summary.tsv"
OUT = ROOT / "genome_5mer_validation.tsv"

LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
CONTEXTS = ["CG", "CHG", "CHH"]

data = defaultdict(list)
seen = set()

with open(SOURCE) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if row["posterior_threshold"] != "0.95":
            continue

        lineage = row["lineage"]
        context = row["context"]
        motif = row["ancestral_5mer"]
        opportunities = int(row["opportunities"])
        ct = int(row["C>T"])

        assert lineage in LINEAGES
        assert context in CONTEXTS
        assert len(motif) == 5
        assert motif[2] == "C"
        assert set(motif) <= set("ACGT")
        assert 0 <= ct <= opportunities

        expected = (
            "CG" if motif[3] == "G"
            else "CHG" if motif[4] == "G"
            else "CHH"
        )
        assert context == expected, (context, motif)

        key = (lineage, context, motif)
        assert key not in seen, key
        seen.add(key)

        rate = ct * 1_000_000 / opportunities if opportunities else 0

        data[(lineage, context)].append({
            "motif": motif,
            "opportunities": opportunities,
            "CT": ct,
            "rate": rate
        })

rows = []

print("===== EXACT 5-MER VALIDATION =====")

for lineage in LINEAGES:
    print("\n" + lineage.upper())

    for context in CONTEXTS:
        motifs = data[(lineage, context)]
        motifs = sorted(motifs, key=lambda x: x["rate"])

        expected_count = {"CG": 64, "CHG": 48, "CHH": 144}[context]
        assert len(motifs) == expected_count, (
            lineage, context, len(motifs)
        )

        lowest = motifs[0]
        highest = motifs[-1]

        print(
            context,
            "motifs:", len(motifs),
            "lowest:", lowest["motif"],
            f'{lowest["rate"]:.2f}',
            "highest:", highest["motif"],
            f'{highest["rate"]:.2f}',
            "fold range:",
            f'{highest["rate"]/lowest["rate"]:.3f}'
            if lowest["rate"] else "NA"
        )

        for m in motifs:
            rows.append({
                "lineage": lineage,
                "context": context,
                "ancestral_5mer": m["motif"],
                "opportunities": m["opportunities"],
                "C>T": m["CT"],
                "CT_rate_per_million": f'{m["rate"]:.6f}'
            })

with open(OUT, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("\nWROTE:", OUT)
print("STATUS: PASS")
