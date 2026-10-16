#!/usr/bin/env python3

import gzip
import csv
from array import array
from collections import defaultdict, Counter

STATE="results/phase1/chr21_7species_ASR.state"
MAP="results/phase1/chr21_7species_alignment_map.tsv.gz"

OUT="results/phase1/chr21_formal_5mer_summary.tsv"
QCOUT="results/phase1/chr21_formal_5mer_QC.tsv"

N=24200593

# Internal nodes from IQ-TREE:
# Node1 = human-Pan ancestor
# Node2 = chimp-bonobo (Pan) ancestor
# Node3 = human-Pan-gorilla ancestor
NODES={
    "Node1":"human_pan",
    "Node2":"pan",
    "Node3":"human_pan_gorilla"
}

code={"A":1,"C":2,"G":3,"T":4}
decode={1:"A",2:"C",3:"G",4:"T"}

qc=Counter()

# Compact storage:
# one byte per ML state + one float per posterior probability.
states={n:bytearray(N) for n in NODES}
probs={n:array('f',[0.0])*N for n in NODES}

print("Reading ASR state file...", flush=True)

with open(STATE,"rt") as fh:

    for line in fh:

        if not line:
            continue

        if line.startswith("#"):
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
        state=f[2]

        qc[f"{node}_rows"]+=1

        # IQ-TREE can occasionally report "-" or another
        # non-ACGT ML state. Treat this as missing/unusable.
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

print("ASR extraction complete.", flush=True)

for node in NODES:
    print(
        node,
        "rows:",qc[f"{node}_rows"],
        "ACGT:",qc[f"{node}_ACGT_calls"],
        "nonACGT:",qc[f"{node}_nonACGT_state"],
        flush=True
    )

comp=str.maketrans("ACGT","TGCA")

def revcomp(s):
    return s.translate(comp)[::-1]

def canonical_mut(anc,derived):

    if anc=="C":
        return "C>"+derived

    if anc=="G":
        # Complement G-centered mutation onto C strand.
        d=derived.translate(comp)
        return "C>"+d

    return None

def context_class(five):

    # five is always oriented so central ancestral base is C
    if len(five)!=5 or five[2]!="C":
        return None

    if five[3]=="G":
        return "CG"

    if five[4]=="G":
        return "CHG"

    return "CHH"


summary=defaultdict(Counter)

lineage_info={
    "human":("Node1","hg38"),
    "chimp":("Node2","panTro5"),
    "bonobo":("Node2","panPan2"),
    "gorilla":("Node3","gorGor5")
}

thresholds=[0.90,0.95,0.99]

print("Streaming coordinate map...", flush=True)

window=[]

with gzip.open(MAP,"rt") as fh:

    r=csv.DictReader(fh,delimiter="\t")

    for row in r:

        site=int(row["site"])-1
        pos=int(row["hg38_pos0"])

        item={
            "site":site,
            "pos":pos,
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

        # Require a genuine consecutive 5-bp hg38 genomic interval.
        if positions != list(range(positions[0],positions[0]+5)):
            qc["noncontiguous_windows"]+=1
            continue

        qc["contiguous_5bp_windows"]+=1

        center=window[2]

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

            # Analyze cytosines on either DNA strand.
            if center_anc not in ("C","G"):
                continue

            five="".join(anc)

            # Canonicalize G-centered sites to C-centered orientation.
            if center_anc=="G":
                five=revcomp(five)

            ctx=context_class(five)

            if ctx is None:
                qc[f"{lineage}_bad_context"]+=1
                continue

            qc[f"{lineage}_ancestral_C_or_G_windows"]+=1

            min_p=min(pp)

            for thr in thresholds:

                # Require every nucleotide of ancestral 5-mer
                # to satisfy posterior threshold.
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


print("Writing summary...", flush=True)

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


print("Writing QC...", flush=True)

with open(QCOUT,"w",newline="") as fh:

    w=csv.writer(fh,delimiter="\t")
    w.writerow(["metric","count"])

    for k in sorted(qc):
        w.writerow([k,qc[k]])


print("DONE", flush=True)
print("Summary:",OUT, flush=True)
print("QC:",QCOUT, flush=True)
