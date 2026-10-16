#!/usr/bin/env python3

import csv
import glob
import itertools
from pathlib import Path

import numpy as np

ROOT = Path("results/phase1_genome")
OUT = Path("results/phase1_strengthening/motif_deviations")
OUT.mkdir(parents=True, exist_ok=True)

LINEAGES = ["human", "chimp", "bonobo", "gorilla"]
CONTEXTS = ["CG", "CHG", "CHH"]
B = 5000
SEED = 20261008
THRESHOLD = 0.95

files = sorted(glob.glob(
    str(ROOT / "chr*_formal_5mer_STRICT_summary.production.tsv")
))

assert len(files) == 22, f"Expected 22 chromosome files, found {len(files)}"

records = {}

for filename in files:
    chrom = Path(filename).name.split("_")[0]

    with open(filename, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if abs(float(row["posterior_threshold"]) - THRESHOLD) > 1e-8:
                continue

            key = (
                chrom,
                row["lineage"],
                row["context"],
                row["ancestral_5mer"]
            )

            if key in records:
                raise ValueError(f"Duplicate row: {key}")

            records[key] = (
                int(row["opportunities"]),
                int(row["C>T"])
            )

chroms = [f"chr{i}" for i in range(1, 23)]
motifs = sorted({
    key[3] for key in records
})

assert len(motifs) == 256, f"Expected 256 motifs, found {len(motifs)}"

motif_context = {}

for key in records:
    motif = key[3]
    context = key[2]
    if motif in motif_context and motif_context[motif] != context:
        raise ValueError(f"Inconsistent context: {motif}")
    motif_context[motif] = context

# Arrays: chromosome x lineage x motif
opp = np.zeros((22, 4, 256), dtype=np.float64)
ct = np.zeros_like(opp)

for c, chrom in enumerate(chroms):
    for l, lineage in enumerate(LINEAGES):
        for m, motif in enumerate(motifs):
            key = (chrom, lineage, motif_context[motif], motif)
            if key not in records:
                raise ValueError(f"Missing record: {key}")
            opp[c, l, m], ct[c, l, m] = records[key]

context_masks = {
    context: np.array(
        [motif_context[m] == context for m in motifs]
    )
    for context in CONTEXTS
}

def relative_effect(opportunities, substitutions):
    """Log2 of motif rate divided by context-wide rate."""
    result = np.full((4, 256), np.nan)

    for context in CONTEXTS:
        mask = context_masks[context]

        base_opp = opportunities[:, mask].sum(axis=1)
        base_ct = substitutions[:, mask].sum(axis=1)

        base_rate = np.divide(
            base_ct, base_opp,
            out=np.full(4, np.nan),
            where=base_opp > 0
        )

        motif_rate = np.divide(
            substitutions[:, mask],
            opportunities[:, mask],
            out=np.full((4, mask.sum()), np.nan),
            where=opportunities[:, mask] > 0
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            result[:, mask] = np.log2(
                motif_rate / base_rate[:, None]
            )

    return result

total_opp = opp.sum(axis=0)
total_ct = ct.sum(axis=0)
observed = relative_effect(total_opp, total_ct)

# Use the SAME chromosome resampling for all lineages,
# preserving cross-lineage chromosome-level dependence.
rng = np.random.default_rng(SEED)
indices = rng.integers(0, 22, size=(B, 22))

bootstrap = np.empty((B, 4, 256), dtype=np.float64)

for b in range(B):
    sampled_opp = opp[indices[b]].sum(axis=0)
    sampled_ct = ct[indices[b]].sum(axis=0)
    bootstrap[b] = relative_effect(sampled_opp, sampled_ct)

def ci(values):
    valid = values[np.isfinite(values)]
    if len(valid) < B * 0.95:
        return np.nan, np.nan
    return tuple(np.quantile(valid, [0.025, 0.975]))

def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

# Table 1: per-lineage motif effects
effects = []

for m, motif in enumerate(motifs):
    for l, lineage in enumerate(LINEAGES):
        lo, hi = ci(bootstrap[:, l, m])

        effects.append({
            "context": motif_context[motif],
            "motif": motif,
            "lineage": lineage,
            "opportunities": int(total_opp[l, m]),
            "CT": int(total_ct[l, m]),
            "relative_log2_effect": observed[l, m],
            "bootstrap_CI_low": lo,
            "bootstrap_CI_high": hi
        })

write_tsv(
    OUT / "lineage_motif_effects.tsv",
    effects,
    list(effects[0])
)

# Table 2: six lineage-pair comparisons
comparisons = []

for l1, l2 in itertools.combinations(range(4), 2):
    for m, motif in enumerate(motifs):
        delta = observed[l1, m] - observed[l2, m]

        boot_delta = (
            bootstrap[:, l1, m] -
            bootstrap[:, l2, m]
        )

        lo, hi = ci(boot_delta)

        comparisons.append({
            "context": motif_context[motif],
            "motif": motif,
            "lineage_1": LINEAGES[l1],
            "lineage_2": LINEAGES[l2],
            "delta_log2_relative_effect": delta,
            "absolute_delta": abs(delta),
            "bootstrap_CI_low": lo,
            "bootstrap_CI_high": hi,
            "CI_excludes_zero": (
                bool(lo > 0 or hi < 0)
                if np.isfinite(lo) and np.isfinite(hi)
                else False
            )
        })

comparisons.sort(
    key=lambda row: (
        -row["absolute_delta"]
        if np.isfinite(row["absolute_delta"])
        else float("inf")
    )
)

write_tsv(
    OUT / "pairwise_motif_deviations.tsv",
    comparisons,
    list(comparisons[0])
)

# Table 3: motifs with largest across-lineage dispersion
ranked = []

for m, motif in enumerate(motifs):
    values = observed[:, m]

    ranked.append({
        "context": motif_context[motif],
        "motif": motif,
        "lineage_effect_range": (
            float(np.ptp(values))
            if np.all(np.isfinite(values))
            else np.nan
        ),
        **{
            f"{lineage}_log2_effect": observed[l, m]
            for l, lineage in enumerate(LINEAGES)
        }
    })

ranked.sort(
    key=lambda row: -row["lineage_effect_range"]
    if np.isfinite(row["lineage_effect_range"])
    else float("inf")
)

write_tsv(
    OUT / "ranked_lineage_deviations.tsv",
    ranked,
    list(ranked[0])
)

# Heatmap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(
    1, 3, figsize=(15, 10),
    gridspec_kw={"width_ratios": [64, 48, 144]}
)

for ax, context in zip(axes, CONTEXTS):
    idx = [
        m for m, motif in enumerate(motifs)
        if motif_context[motif] == context
    ]

    matrix = observed[:, idx]
    matrix = matrix - np.nanmean(
        matrix, axis=0, keepdims=True
    )

    image = ax.imshow(
        matrix,
        aspect="auto",
        cmap="coolwarm",
        vmin=-0.5,
        vmax=0.5,
        interpolation="nearest"
    )

    ax.set_title(context)
    ax.set_yticks(range(4))
    ax.set_yticklabels(LINEAGES)
    ax.set_xticks([])

fig.colorbar(
    image, ax=axes.tolist(),
    shrink=0.55,
    label="Lineage deviation from motif mean (log2)"
)

fig.suptitle(
    "Lineage-specific deviations in ancestral 5-mer C>T patterns"
)

fig.subplots_adjust(
    left=0.09, right=0.87,
    bottom=0.08, top=0.90, wspace=0.22
)

fig.savefig(
    OUT / "lineage_motif_deviation_heatmap.png",
    dpi=300
)
fig.savefig(
    OUT / "lineage_motif_deviation_heatmap.pdf"
)
plt.close(fig)

print("ANALYSIS COMPLETE")
print("Chromosomes:", len(chroms))
print("Motifs:", len(motifs))
print("Bootstrap replicates:", B)

print("\nTOP 15 LINEAGE-PAIR DEVIATIONS")
for row in comparisons[:15]:
    print(
        row["context"],
        row["motif"],
        row["lineage_1"],
        row["lineage_2"],
        f"delta={row['delta_log2_relative_effect']:.4f}",
        f"CI=[{row['bootstrap_CI_low']:.4f},"
        f"{row['bootstrap_CI_high']:.4f}]",
        f"excludes_zero={row['CI_excludes_zero']}"
    )

print("\nOutput:", OUT)
