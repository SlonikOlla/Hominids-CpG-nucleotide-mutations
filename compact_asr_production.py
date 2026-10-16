#!/usr/bin/env python3

import gzip
import csv
import sys
import os
from array import array
from collections import defaultdict, Counter

if len(sys.argv) != 2:
    sys.exit("Usage: compact_asr_production.py chr22")

CHR=sys.argv[1]

ROOT="/scratch/vasilisa/hominid_methyl_evolution"
WORK=f"{ROOT}/work/{CHR}"

FINAL="/project/def-vasilisa/hominid_methyl_evolution/results/phase1_genome"
os.makedirs(FINAL,exist_ok=True)

STATE=f"{WORK}/{CHR}_7species_ASR.state"
MAP=f"{WORK}/{CHR}_7species_production_map.tsv.gz"
PHY=f"{WORK}/{CHR}_7species.phy"

OUT=f"{FINAL}/{CHR}_formal_5mer_STRICT_summary.production.tsv"
QCOUT=f"{FINAL}/{CHR}_formal_5mer_STRICT_QC.production.tsv"
SITEOUT=f"{FINAL}/{CHR}_ancestral_sites_ge0.95.tsv.gz"

# ------------------------------------------------------------
# Alignment size
# ------------------------------------------------------------

with open(PHY) as fh:
    header=fh.readline().split()

if len(header)<2:
    sys.exit(f"ERROR: malformed PHYLIP header in {PHY}")

N=int(header[1])

print("Chromosome:",CHR,flush=True)
print("Alignment sites:",N,flush=True)

# ------------------------------------------------------------
# Internal nodes
# ------------------------------------------------------------

# Node1 = human-Pan ancestor
# Node2 = chimp-bonobo ancestor
# Node3 = human-Pan-gorilla ancestor

NODES={
    "Node1":"human_pan",
    "Node2":"pan",
    "Node3":"human_pan_gorilla"
}

code={"A":1,"C":2,"G":3,"T":4}
decode={1:"A",2:"C",3:"G",4:"T"}

qc=Counter()

# One byte per ML state + one float per posterior probability.
states={n:bytearray(N) for n in NODES}
probs={n:array('f',[0.0])*N for n in NODES}

# ------------------------------------------------------------
# Read ASR
# ------------------------------------------------------------

print("Reading ASR state file...",flush=True)

with open(STATE,"rt") as fh:

    for line in fh:

        if not line or line.startswith("#"):
            continue

        if line.startswith("Node\t") or line.startswith("Node "):
            continue

        f=line.rstrip().split()

        if len(f)<7:
            qc["ASR_malformed_rows"]+=1
            continue

        node=f[0]

        if node not in NODES:
            continue

        site=int(f[1])-1

        if site < 0 or site >= N:
            qc["ASR_site_out_of_range"]+=1
            continue

        state=f[2]

        qc[f"{node}_rows"]+=1

        if state not in code:
            qc[f"{node}_nonACGT_state"]+=1
            continue

        p={
            "A":float(f[3]),
            "C":float(f[4]),
            "G":float(f[5]),
            "T":float(f[6])
        }

        states[node][site]=code[state]
        probs[node][site]=p[state]

        qc[f"{node}_ACGT_calls"]+=1

print("ASR extraction complete.",flush=True)

for node in NODES:
    print(
        node,
        "rows:",qc[f"{node}_rows"],
        "ACGT:",qc[f"{node}_ACGT_calls"],
        "nonACGT:",qc[f"{node}_nonACGT_state"],
        flush=True
    )

# ------------------------------------------------------------
# Sequence helpers
# ------------------------------------------------------------

comp=str.maketrans("ACGT","TGCA")

def revcomp(s):
    return s.translate(comp)[::-1]

def canonical_mut(anc,derived):

    if anc=="C":
        if derived=="C":
            return "."
        return "C>"+derived

    if anc=="G":
        if derived=="G":
            return "."
        return "C>"+derived.translate(comp)

    return None

def canonical_base(anc,derived):

    if anc=="C":
        return derived

    if anc=="G":
        return derived.translate(comp)

    return None

def context_class(five):

    # five must already be C-centered.
    if len(five)!=5 or five[2]!="C":
        return None

    if five[3]=="G":
        return "CG"

    if five[4]=="G":
        return "CHG"

    return "CHH"

# ------------------------------------------------------------
# Analysis definitions
# ------------------------------------------------------------

summary=defaultdict(Counter)

lineage_info={
    "human":("Node1","hg38"),
    "chimp":("Node2","panTro5"),
    "bonobo":("Node2","panPan2"),
    "gorilla":("Node3","gorGor5")
}

thresholds=[0.90,0.95,0.99]
PRIMARY=0.95

# ------------------------------------------------------------
# Site-level output
#
# One row per hg38 genomic coordinate/window.
# A lineage's fields are "." if that lineage does not satisfy
# the primary >=0.95 complete ancestral-5mer criterion.
# ------------------------------------------------------------

site_fields=[
    "chr",
    "hg38_pos0",
    "hg38_pos1",
    "site",
    "maf_block",
    "aln_col"
]

for lineage in lineage_info:
    site_fields.extend([
        f"{lineage}_ancestral_base",
        f"{lineage}_center_PP",
        f"{lineage}_ancestral_5mer",
        f"{lineage}_min_5mer_PP",
        f"{lineage}_context",
        f"{lineage}_descendant_base",
        f"{lineage}_descendant_base_canonical",
        f"{lineage}_mutation"
    ])

print("Streaming coordinate map...",flush=True)

window=[]

with gzip.open(
    SITEOUT,"wt",compresslevel=6,newline=""
) as site_fh, gzip.open(MAP,"rt") as map_fh:

    site_writer=csv.DictWriter(
        site_fh,
        fieldnames=site_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    site_writer.writeheader()

    r=csv.DictReader(map_fh,delimiter="\t")

    map_rows=0

    for row in r:

        map_rows+=1

        site=int(row["site"])-1
        pos=int(row["hg38_pos0"])

        if site < 0 or site >= N:
            qc["map_site_out_of_range"]+=1
            continue

        item={
            "site":site,
            "pos":pos,
            "block":int(row["maf_block"]),
            "aln_col":int(row["aln_col"]),
            "bases":{
                "hg38":row["hg38"],
                "panTro5":row["panTro5"],
                "panPan2":row["panPan2"],
                "gorGor5":row["gorGor5"]
            }
        }

        window.append(item)

        if len(window)<5:
            continue

        if len(window)>5:
            window.pop(0)

        positions=[x["pos"] for x in window]
        blocks=[x["block"] for x in window]
        cols=[x["aln_col"] for x in window]

        # Strict production definition.
        if len(set(blocks)) != 1:
            qc["window_crosses_MAF_block"]+=1
            continue

        if cols != list(range(cols[0],cols[0]+5)):
            qc["window_skips_alignment_column"]+=1
            continue

        if positions != list(range(positions[0],positions[0]+5)):
            qc["window_nonconsecutive_hg38"]+=1
            continue

        qc["strict_contiguous_5bp_windows"]+=1

        center=window[2]

        # Site-level row assembled across lineages.
        site_record={
            "chr":CHR,
            "hg38_pos0":center["pos"],
            "hg38_pos1":center["pos"]+1,
            "site":center["site"]+1,
            "maf_block":center["block"],
            "aln_col":center["aln_col"]
        }

        any_primary=False

        for lineage in lineage_info:
            for suffix in [
                "ancestral_base",
                "center_PP",
                "ancestral_5mer",
                "min_5mer_PP",
                "context",
                "descendant_base",
                "descendant_base_canonical",
                "mutation"
            ]:
                site_record[f"{lineage}_{suffix}"]="."

        for lineage,(node,focal) in lineage_info.items():

            anc=[]
            pp=[]
            valid=True

            for x in window:

                s=x["site"]
                c=states[node][s]

                if c==0:
                    valid=False
                    break

                anc.append(decode[c])
                pp.append(probs[node][s])

            if not valid:
                qc[f"{lineage}_window_missing_ASR"]+=1
                continue

            focal_base=center["bases"][focal]

            if focal_base not in "ACGT":
                qc[f"{lineage}_bad_focal"]+=1
                continue

            center_anc=anc[2]

            if center_anc not in ("C","G"):
                continue

            five="".join(anc)

            if center_anc=="G":
                five=revcomp(five)

            ctx=context_class(five)

            if ctx is None:
                qc[f"{lineage}_bad_context"]+=1
                continue

            qc[f"{lineage}_ancestral_C_or_G_windows"]+=1

            min_p=min(pp)
            center_p=pp[2]

            # Validated summary calculations remain unchanged.
            for thr in thresholds:

                if min_p < thr:
                    continue

                key=(lineage,thr,ctx,five)

                summary[key]["opportunities"]+=1
                qc[f"{lineage}_windows_ge_{thr:.2f}"]+=1

                if focal_base != center_anc:

                    mut=canonical_mut(center_anc,focal_base)

                    if mut in ("C>T","C>A","C>G"):

                        summary[key][mut]+=1
                        summary[key]["total"]+=1
                        qc[f"{lineage}_substitutions_ge_{thr:.2f}"]+=1

            # Primary-threshold site table.
            if min_p >= PRIMARY:

                any_primary=True

                mut=canonical_mut(center_anc,focal_base)
                canonical_desc=canonical_base(center_anc,focal_base)

                site_record[f"{lineage}_ancestral_base"]="C"
                site_record[f"{lineage}_center_PP"]=f"{center_p:.6f}"
                site_record[f"{lineage}_ancestral_5mer"]=five
                site_record[f"{lineage}_min_5mer_PP"]=f"{min_p:.6f}"
                site_record[f"{lineage}_context"]=ctx
                site_record[f"{lineage}_descendant_base"]=focal_base
                site_record[f"{lineage}_descendant_base_canonical"]=canonical_desc
                site_record[f"{lineage}_mutation"]=mut

        # Write once if >=1 lineage provides a primary opportunity.
        if any_primary:
            site_writer.writerow(site_record)
            qc["coordinate_rows_ge_0.95_any_lineage"]+=1

    qc["map_rows"]=map_rows

# ------------------------------------------------------------
# Internal consistency checks
# ------------------------------------------------------------

if qc["map_rows"] != N:
    raise RuntimeError(
        f"Map/PHYLIP mismatch: map={qc['map_rows']} PHYLIP={N}"
    )

for node in NODES:
    if qc[f"{node}_rows"] != N:
        raise RuntimeError(
            f"{node} row-count mismatch: "
            f"{qc[f'{node}_rows']} versus {N}"
        )

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("Writing summary...",flush=True)

with open(OUT,"w",newline="") as fh:

    w=csv.writer(fh,delimiter="\t")

    w.writerow([
        "lineage",
        "posterior_threshold",
        "context",
        "ancestral_5mer",
        "opportunities",
        "C>T",
        "C>A",
        "C>G",
        "total_substitutions",
        "rate_per_million",
        "CT_rate_per_million"
    ])

    for key in sorted(summary):

        lineage,thr,ctx,five=key
        x=summary[key]

        opp=x["opportunities"]
        total=x["total"]

        rate=(total/opp*1e6) if opp else 0
        ctr=(x["C>T"]/opp*1e6) if opp else 0

        w.writerow([
            lineage,
            f"{thr:.2f}",
            ctx,
            five,
            opp,
            x["C>T"],
            x["C>A"],
            x["C>G"],
            total,
            f"{rate:.6f}",
            f"{ctr:.6f}"
        ])

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

print("Writing QC...",flush=True)

with open(QCOUT,"w",newline="") as fh:

    w=csv.writer(fh,delimiter="\t")
    w.writerow(["metric","count"])

    for k in sorted(qc):
        w.writerow([k,qc[k]])

print("DONE",flush=True)
print("Summary:",OUT,flush=True)
print("QC:",QCOUT,flush=True)
print("Site table:",SITEOUT,flush=True)
