#!/usr/bin/env python3

import gzip

MAF="data/maf/chr21.maf.gz"
OUT="results/phase1/chr21_7species.phy"

species=[
    "hg38",
    "panTro5",
    "panPan2",
    "gorGor5",
    "ponAbe2",
    "nomLeu3",
    "rheMac8"
]

VALID=set("ACGT")

seqs={s:[] for s in species}

def sp(src):
    return src.split(".",1)[0]

def process(block):
    if not all(s in block for s in species):
        return 0

    lengths={len(block[s]) for s in species}
    if len(lengths)!=1:
        return 0

    L=lengths.pop()
    kept=0

    for i in range(L):
        col=[block[s][i] for s in species]

        if all(x in VALID for x in col):
            for s,x in zip(species,col):
                seqs[s].append(x)
            kept+=1

    return kept

block={}
kept_total=0
blocks=0
usable_blocks=0

with gzip.open(MAF,"rt") as fh:

    for line in fh:

        if line.startswith("a "):

            if block:
                blocks+=1
                n=process(block)
                kept_total+=n
                if n:
                    usable_blocks+=1

            block={}

        elif line.startswith("s "):

            f=line.rstrip().split()
            s=sp(f[1])

            if s in species:
                block[s]=f[6].upper()

    if block:
        blocks+=1
        n=process(block)
        kept_total+=n
        if n:
            usable_blocks+=1

print("Blocks:",blocks)
print("Usable blocks:",usable_blocks)
print("Complete 7-species columns:",kept_total)

with open(OUT,"w") as out:

    out.write(f"{len(species)} {kept_total}\n")

    for s in species:
        sequence="".join(seqs[s])
        assert len(sequence)==kept_total
        out.write(f"{s:<12} {sequence}\n")

print("Wrote:",OUT)
