#!/usr/bin/env python3

import gzip
import csv
from collections import Counter, defaultdict

MAF = "data/maf/chr21.maf.gz"
OUT_SUM = "results/phase1/chr21_human_context_summary.tsv"
OUT_SITES = "results/phase1/chr21_human_derived_sites.tsv.gz"
OUT_QC = "results/phase1/chr21_human_pilot_QC.tsv"

TARGETS = [
    "hg38",
    "panTro5",
    "panPan2",
    "gorGor5",
    "ponAbe2",
    "nomLeu3",
    "rheMac8",
]

VALID = set("ACGT")

qc = Counter()
summary = defaultdict(Counter)

def species(src):
    return src.split(".", 1)[0]

def parse_blocks(path):
    block = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("a "):
                if block:
                    yield block
                block = {}
            elif line.startswith("s "):
                f = line.rstrip().split()
                sp = species(f[1])
                if sp in TARGETS:
                    block[sp] = {
                        "src": f[1],
                        "start": int(f[2]),
                        "size": int(f[3]),
                        "strand": f[4],
                        "srcSize": int(f[5]),
                        "seq": f[6].upper(),
                    }
        if block:
            yield block

def context_class(tri):
    # tri is oriented so focal ancestral C is tri[0]
    if len(tri) < 2:
        return None
    if tri[0] != "C":
        return None
    if tri[1] == "G":
        return "CG"
    if len(tri) >= 3 and tri[2] == "G":
        return "CHG"
    return "CHH"

def human_positions(rec):
    """
    Return hg38 coordinate for every alignment column.
    MAF start is zero-based. hg38 should be '+' in this hg38-centered MAF,
    but handle both strands defensively.
    """
    seq = rec["seq"]
    pos = rec["start"]
    coords = []
    if rec["strand"] == "+":
        p = pos
        for b in seq:
            if b == "-":
                coords.append(None)
            else:
                coords.append(p)
                p += 1
    else:
        # Convert reverse-strand MAF coordinates to forward coordinates.
        p = rec["srcSize"] - rec["start"] - 1
        for b in seq:
            if b == "-":
                coords.append(None)
            else:
                coords.append(p)
                p -= 1
    return coords

with gzip.open(OUT_SITES, "wt", newline="") as sf:
    writer = csv.writer(sf, delimiter="\t")
    writer.writerow([
        "chr", "hg38_pos_1based",
        "ancestral_context", "context_class",
        "human_base", "ancestral_base",
        "substitution",
        "chimp", "bonobo", "gorilla",
        "orangutan", "gibbon", "rhesus"
    ])

    for block_no, block in enumerate(parse_blocks(MAF), 1):
        qc["blocks_total"] += 1

        # Require human, chimp and bonobo for this conservative pilot.
        if not all(x in block for x in ("hg38", "panTro5", "panPan2")):
            qc["blocks_missing_core"] += 1
            continue

        lengths = {len(x["seq"]) for x in block.values()}
        if len(lengths) != 1:
            qc["blocks_bad_alignment_length"] += 1
            continue

        L = lengths.pop()
        coords = human_positions(block["hg38"])

        # Build an inferred ancestral base for every alignment column.
        anc = [None] * L

        for i in range(L):
            h = block["hg38"]["seq"][i]
            c = block["panTro5"]["seq"][i]
            b = block["panPan2"]["seq"][i]

            if h not in VALID or c not in VALID or b not in VALID:
                qc["columns_core_gap_or_N"] += 1
                continue

            # Chimp + bonobo must agree.
            if c != b:
                qc["columns_chimp_bonobo_disagree"] += 1
                continue

            support = 0
            for sp in ("gorGor5", "ponAbe2", "nomLeu3"):
                if sp in block:
                    x = block[sp]["seq"][i]
                    if x == c:
                        support += 1

            # At least one additional hominoid must support Pan state.
            if support < 1:
                qc["columns_no_extra_hominoid_support"] += 1
                continue

            anc[i] = c
            qc["columns_ancestral_callable"] += 1

        # Examine ancestral C sites and derive context using the next
        # two callable ancestral nucleotides in ungapped ancestral sequence.
        callable_cols = [i for i, x in enumerate(anc) if x in VALID]

        for k, i in enumerate(callable_cols):
            if anc[i] != "C":
                continue

            qc["ancestral_C_sites"] += 1

            # Need next two contiguous alignment columns to avoid crossing
            # gaps/uncertain ancestral positions.
            if i + 2 >= L:
                qc["C_no_context"] += 1
                continue

            if anc[i+1] not in VALID or anc[i+2] not in VALID:
                qc["C_no_context"] += 1
                continue

            # Require true adjacency in human coordinate as an additional
            # conservative guard against indels.
            if coords[i] is None or coords[i+1] is None or coords[i+2] is None:
                qc["C_noncontiguous_human"] += 1
                continue

            if not (coords[i+1] == coords[i] + 1 and
                    coords[i+2] == coords[i] + 2):
                qc["C_noncontiguous_human"] += 1
                continue

            tri = anc[i] + anc[i+1] + anc[i+2]
            cls = context_class(tri)

            if cls is None:
                qc["C_bad_context"] += 1
                continue

            hb = block["hg38"]["seq"][i]
            if hb not in VALID:
                qc["human_uncallable_at_C"] += 1
                continue

            summary[cls]["opportunities"] += 1
            qc["C_context_callable"] += 1

            if hb != "C":
                sub = f"C>{hb}"
                summary[cls]["substitutions"] += 1
                summary[cls][sub] += 1
                qc["human_derived_C_substitutions"] += 1

                vals = {}
                for sp, label in [
                    ("panTro5","chimp"),
                    ("panPan2","bonobo"),
                    ("gorGor5","gorilla"),
                    ("ponAbe2","orangutan"),
                    ("nomLeu3","gibbon"),
                    ("rheMac8","rhesus"),
                ]:
                    vals[label] = (
                        block[sp]["seq"][i]
                        if sp in block else "."
                    )

                writer.writerow([
                    block["hg38"]["src"].split(".",1)[1],
                    coords[i] + 1,
                    tri,
                    cls,
                    hb,
                    "C",
                    sub,
                    vals["chimp"],
                    vals["bonobo"],
                    vals["gorilla"],
                    vals["orangutan"],
                    vals["gibbon"],
                    vals["rhesus"],
                ])

with open(OUT_SUM, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow([
        "lineage", "context",
        "opportunities",
        "C>T", "C>A", "C>G",
        "total_substitutions",
        "rate_per_million"
    ])

    for cls in ("CG", "CHG", "CHH"):
        x = summary[cls]
        opp = x["opportunities"]
        total = x["substitutions"]
        rate = total / opp * 1_000_000 if opp else 0
        writer.writerow([
            "human", cls, opp,
            x["C>T"], x["C>A"], x["C>G"],
            total,
            f"{rate:.6f}"
        ])

with open(OUT_QC, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["metric", "count"])
    for key in sorted(qc):
        writer.writerow([key, qc[key]])

print("DONE")
print("Summary:", OUT_SUM)
print("Sites:", OUT_SITES)
print("QC:", OUT_QC)
