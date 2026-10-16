#!/usr/bin/env python3

import gzip
import csv

MAF="data/maf/chr21.maf.gz"
OUT="results/phase1/chr21_7species_alignment_map.tsv.gz"

species=[
    "hg38","panTro5","panPan2","gorGor5",
    "ponAbe2","nomLeu3","rheMac8"
]

VALID=set("ACGT")

def sp(src):
    return src.split(".",1)[0]

site=0
block={}

with gzip.open(MAF,"rt") as fh, \
     gzip.open(OUT,"wt",newline="") as out:

    w=csv.writer(out,delimiter="\t")

    w.writerow([
        "site","chr","hg38_pos0",
        *species
    ])

    def process(block):
        global site

        if not all(s in block for s in species):
            return

        lengths={len(block[s]["seq"]) for s in species}
        if len(lengths)!=1:
            return

        L=lengths.pop()

        h=block["hg38"]

        # hg38 is expected on + strand in hg38-centered MAF.
        # Track coordinate through alignment columns.
        pos=h["start"]

        for i in range(L):

            hb=h["seq"][i]

            col=[block[s]["seq"][i] for s in species]

            if hb!="-":
                current_pos=pos
                pos+=1
            else:
                current_pos=None

            if current_pos is None:
                continue

            if not all(x in VALID for x in col):
                continue

            site+=1

            w.writerow([
                site,
                h["src"].split(".",1)[1],
                current_pos,
                *col
            ])

    for line in fh:

        if line.startswith("a "):

            if block:
                process(block)

            block={}

        elif line.startswith("s "):

            f=line.rstrip().split()

            s=sp(f[1])

            if s in species:

                block[s]={
                    "src":f[1],
                    "start":int(f[2]),
                    "size":int(f[3]),
                    "strand":f[4],
                    "srcSize":int(f[5]),
                    "seq":f[6].upper()
                }

    if block:
        process(block)

print("Mapped sites:",site)
print("Output:",OUT)
