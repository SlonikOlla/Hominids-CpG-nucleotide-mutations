#!/usr/bin/env python3

import gzip
import csv
import sys
import os

if len(sys.argv) != 2:
    sys.exit("Usage: make_production_map.py chr22")

CHR=sys.argv[1]

ROOT="/scratch/vasilisa/hominid_methyl_evolution"
MAF=f"{ROOT}/maf/{CHR}.maf.gz"
WORK=f"{ROOT}/work/{CHR}"
OUT=f"{WORK}/{CHR}_7species_production_map.tsv.gz"

os.makedirs(WORK,exist_ok=True)

SPECIES=[
    "hg38",
    "panTro5",
    "panPan2",
    "gorGor5",
    "ponAbe2",
    "nomLeu3",
    "rheMac8"
]

def species_name(src):
    return src.split(".",1)[0]

site=0
block_id=0
usable_blocks=0

with gzip.open(MAF,"rt") as fh, \
     gzip.open(OUT,"wt",newline="") as out:

    w=csv.writer(out,delimiter="\t")

    w.writerow([
        "site","chr","hg38_pos0",
        "maf_block","aln_col",*SPECIES
    ])

    block=[]

    def process_block(block,block_id):
        global site,usable_blocks

        seqs={}
        meta={}

        for line in block:
            if not line.startswith("s "):
                continue

            f=line.split()

            if len(f)<7:
                continue

            src=f[1]
            sp=species_name(src)

            if sp not in SPECIES or sp in seqs:
                continue

            start=int(f[2])
            size=int(f[3])
            strand=f[4]
            srcSize=int(f[5])
            seq=f[6].upper()

            seqs[sp]=seq
            meta[sp]=(src,start,size,strand,srcSize)

        if any(sp not in seqs for sp in SPECIES):
            return

        lengths={len(seqs[sp]) for sp in SPECIES}

        if len(lengths)!=1:
            raise RuntimeError(
                f"Alignment length mismatch in block {block_id}"
            )

        hg_src,hg_start,hg_size,hg_strand,hg_srcSize=meta["hg38"]

        if hg_strand=="+":
            hg_pos=hg_start
        else:
            hg_pos=hg_srcSize-hg_start-1

        block_used=False
        L=len(seqs["hg38"])

        for col in range(L):
            hgbase=seqs["hg38"][col]

            current_pos=None
            if hgbase!="-":
                current_pos=hg_pos

            bases=[seqs[sp][col] for sp in SPECIES]

            if all(b in "ACGT" for b in bases):
                if current_pos is None:
                    raise RuntimeError(
                        f"Complete column without hg38 coordinate: "
                        f"block={block_id} col={col}"
                    )

                site+=1
                block_used=True

                w.writerow([
                    site,
                    hg_src.split(".",1)[1]
                    if "." in hg_src else hg_src,
                    current_pos,
                    block_id,
                    col,
                    *bases
                ])

            if hgbase!="-":
                if hg_strand=="+":
                    hg_pos+=1
                else:
                    hg_pos-=1

        if block_used:
            usable_blocks+=1

    for line in fh:
        if line.startswith("a "):
            if block:
                process_block(block,block_id)

            block_id+=1
            block=[line]

        elif block:
            block.append(line)

    if block:
        process_block(block,block_id)

print("Chromosome:",CHR)
print("MAF blocks:",block_id)
print("Usable blocks:",usable_blocks)
print("Mapped sites:",site)
print("Output:",OUT)
