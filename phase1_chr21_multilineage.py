#!/usr/bin/env python3

import gzip
import csv
from collections import Counter, defaultdict

MAF="data/maf/chr21.maf.gz"
OUT="results/phase1/chr21_multilineage_summary.tsv"
QCOUT="results/phase1/chr21_multilineage_QC.tsv"

TARGETS=[
    "hg38","panTro5","panPan2","gorGor5",
    "ponAbe2","nomLeu3","rheMac8"
]
VALID=set("ACGT")

# For each focal lineage:
# infer ancestral state from phylogenetically appropriate relatives.
#
# Rules intentionally conservative for pilot.
RULES = {
    "human": {
        "focal":"hg38",
        "pair":["panTro5","panPan2"],
        "support":["gorGor5","ponAbe2","nomLeu3"]
    },
    "chimp": {
        "focal":"panTro5",
        "pair":["panPan2","hg38"],
        "support":["gorGor5","ponAbe2","nomLeu3"]
    },
    "bonobo": {
        "focal":"panPan2",
        "pair":["panTro5","hg38"],
        "support":["gorGor5","ponAbe2","nomLeu3"]
    },
    "gorilla": {
        "focal":"gorGor5",
        "pair":["hg38","panTro5"],
        "support":["panPan2","ponAbe2","nomLeu3"]
    }
}

def sp(src):
    return src.split(".",1)[0]

def blocks(path):
    block={}
    with gzip.open(path,"rt") as fh:
        for line in fh:
            if line.startswith("a "):
                if block:
                    yield block
                block={}
            elif line.startswith("s "):
                f=line.rstrip().split()
                s=sp(f[1])
                if s in TARGETS:
                    block[s]=f[6].upper()
        if block:
            yield block

def classify(tri):
    if tri[0]!="C":
        return None
    if tri[1]=="G":
        return "CG"
    if tri[2]=="G":
        return "CHG"
    return "CHH"

summary=defaultdict(Counter)
qc=Counter()

for block in blocks(MAF):
    qc["blocks"]+=1

    lengths={len(x) for x in block.values()}
    if len(lengths)!=1:
        qc["bad_length"]+=1
        continue

    L=lengths.pop()

    for lineage,rule in RULES.items():

        focal=rule["focal"]
        p1,p2=rule["pair"]

        if focal not in block or p1 not in block or p2 not in block:
            qc[f"{lineage}_missing_core"]+=1
            continue

        anc=[None]*L

        for i in range(L):

            fb=block[focal][i]
            a=block[p1][i]
            b=block[p2][i]

            if fb not in VALID or a not in VALID or b not in VALID:
                continue

            # Primary relatives must agree.
            if a!=b:
                continue

            # At least one deeper hominoid supports state.
            support=0
            for s in rule["support"]:
                if s in block and block[s][i]==a:
                    support+=1

            if support<1:
                continue

            anc[i]=a

        for i in range(L-2):

            if anc[i]!="C":
                continue

            # Strict contiguous ancestral context.
            if anc[i+1] not in VALID or anc[i+2] not in VALID:
                continue

            tri=anc[i]+anc[i+1]+anc[i+2]
            ctx=classify(tri)

            if ctx is None:
                continue

            fb=block[focal][i]

            if fb not in VALID:
                continue

            summary[(lineage,ctx)]["opportunities"]+=1

            if fb!="C":
                summary[(lineage,ctx)]["total"]+=1
                summary[(lineage,ctx)][f"C>{fb}"]+=1

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow([
        "lineage","context","opportunities",
        "C>T","C>A","C>G","total",
        "total_rate_per_million",
        "CT_rate_per_million"
    ])

    for lineage in ["human","chimp","bonobo","gorilla"]:
        for ctx in ["CG","CHG","CHH"]:

            x=summary[(lineage,ctx)]
            opp=x["opportunities"]
            total=x["total"]

            tr=total/opp*1e6 if opp else 0
            ctr=x["C>T"]/opp*1e6 if opp else 0

            w.writerow([
                lineage,ctx,opp,
                x["C>T"],x["C>A"],x["C>G"],total,
                f"{tr:.4f}",
                f"{ctr:.4f}"
            ])

with open(QCOUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["metric","count"])
    for k in sorted(qc):
        w.writerow([k,qc[k]])

print("DONE")
print(OUT)
