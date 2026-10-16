#!/usr/bin/env python3

import gzip
import csv

MAF="data/maf/chr21.maf.gz"
OUT="results/phase1/chr21_7species_production_map.tsv.gz"

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
        "site",
        "chr",
        "hg38_pos0",
        "maf_block",
        "aln_col",
        *SPECIES
    ])

    block=[]

    def process_block(block, block_id):
        global site, usable_blocks

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

            if sp not in SPECIES:
                continue

            # Keep first sequence for a species in the block.
            if sp in seqs:
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

        # hg38-centered chromosome MAFs should normally have hg38 on +.
        # Support either strand explicitly.
        hg_pos=None

        if hg_strand=="+":
            hg_pos=hg_start
        else:
            hg_pos=hg_srcSize-hg_start-1

        block_used=False

        L=len(seqs["hg38"])

        for col in range(L):

            hgbase=seqs["hg38"][col]

            # Record current genomic coordinate BEFORE advancing it.
            current_pos=None

            if hgbase!="-":
                current_pos=hg_pos

            bases=[seqs[sp][col] for sp in SPECIES]

            # Exact same complete-case criterion as the original
            # 7-species alignment/map.
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
                    hg_src.split(".",1)[1] if "." in hg_src else hg_src,
                    current_pos,
                    block_id,
                    col,
                    *bases
                ])

            # Advance hg38 genomic coordinate only on non-gap columns.
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


print("MAF blocks:",block_id)
print("Usable blocks:",usable_blocks)
print("Mapped sites:",site)
print("Output:",OUT)

if site != 24200593:
    raise RuntimeError(
        f"Expected 24200593 sites, obtained {site}"
    )
