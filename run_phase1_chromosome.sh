#!/bin/bash
#SBATCH --job-name=phase1_chr
#SBATCH --account=def-vasilisa
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/vasilisa/hominid_methyl_evolution/logs/phase1_%x_%j.out
#SBATCH --error=/scratch/vasilisa/hominid_methyl_evolution/logs/phase1_%x_%j.err

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 chr20"
    exit 1
fi

CHR="$1"

PROJECT=/project/def-vasilisa/hominid_methyl_evolution
SCRATCH=/scratch/vasilisa/hominid_methyl_evolution
WORK="$SCRATCH/work/$CHR"
FINAL="$PROJECT/results/phase1_genome"

cd "$PROJECT"

mkdir -p "$WORK" "$FINAL" "$SCRATCH/logs"

echo "=================================================="
echo "PHASE 1 PRODUCTION"
echo "Chromosome: $CHR"
echo "Started: $(date)"
echo "=================================================="

echo
echo "===== 1. ALIGNMENT ====="

python3 scripts/make_7species_alignment.py "$CHR"

echo
echo "===== 2. PRODUCTION MAP ====="

python3 scripts/make_production_map.py "$CHR"

PHY="$WORK/${CHR}_7species.phy"
MAP="$WORK/${CHR}_7species_production_map.tsv.gz"

echo
echo "===== 3. ALIGNMENT/MAP CROSS-CHECK ====="

PHY_N=$(head -1 "$PHY" | awk '{print $2}')
MAP_N=$(gzip -cd "$MAP" | awk 'END{print NR-1}')

echo "PHYLIP sites: $PHY_N"
echo "Map sites:    $MAP_N"

if [[ "$PHY_N" != "$MAP_N" ]]; then
    echo "ERROR: PHYLIP/map mismatch"
    exit 1
fi

echo "PASS: PHYLIP/map site counts identical"

echo
echo "===== 4. IQ-TREE ASR ====="

module load iq-tree/2.3.6

cp "$PROJECT/data/metadata/hominid_7species.tree" \
   "$WORK/hominid_7species.tree"

cd "$WORK"

iqtree2 \
    -s "${CHR}_7species.phy" \
    -m GTR+F+G4 \
    -te hominid_7species.tree \
    -asr \
    -pre "${CHR}_7species_ASR" \
    -nt 1 \
    -redo

STATE="$WORK/${CHR}_7species_ASR.state"

if [[ ! -s "$STATE" ]]; then
    echo "ERROR: ASR state file missing"
    exit 1
fi

echo
echo "===== 5. VERIFY NODE TOPOLOGY ====="

TREEFILE="$WORK/${CHR}_7species_ASR.treefile"

cat "$TREEFILE"

python3 - "$TREEFILE" <<'PY'
import sys,re

tree=open(sys.argv[1]).read().strip().replace("\n","")

tests={
    "Node2":r"\(panTro5:[^,()]+,panPan2:[^,()]+\)Node2:",
    "Node1":r"\(hg38:[^,()]+,\(panTro5:[^,()]+,panPan2:[^,()]+\)Node2:[^,()]+,",
    "Node3":r"\(gorGor5:[^,()]+,\(ponAbe2:"
}

for node,pat in tests.items():
    if not re.search(pat,tree):
        raise SystemExit(f"ERROR: expected topology for {node} not found")

print("PASS: expected Node1/Node2/Node3 topology detected")
PY

echo
echo "===== 6. REDUCE ASR ====="

cd "$PROJECT"

python3 scripts/compact_asr_production.py "$CHR"

SUMMARY="$FINAL/${CHR}_formal_5mer_STRICT_summary.production.tsv"
QC="$FINAL/${CHR}_formal_5mer_STRICT_QC.production.tsv"
SITES="$FINAL/${CHR}_ancestral_sites_ge0.95.tsv.gz"

echo
echo "===== 7. FINAL OUTPUT CHECK ====="

for f in "$SUMMARY" "$QC" "$SITES"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing/empty final output: $f"
        exit 1
    fi
done

python3 - "$QC" "$PHY_N" <<'PY'
import csv,sys

qcfile=sys.argv[1]
N=int(sys.argv[2])

with open(qcfile) as f:
    q={
        x["metric"]:int(x["count"])
        for x in csv.DictReader(f,delimiter="\t")
    }

required=[
    "map_rows",
    "strict_contiguous_5bp_windows",
    "coordinate_rows_ge_0.95_any_lineage"
]

for k in required:
    if k not in q:
        raise SystemExit(f"ERROR: required QC metric missing: {k}")

if q["map_rows"] != N:
    raise SystemExit(
        f"ERROR: QC map_rows={q['map_rows']} but alignment N={N}"
    )

for node in ["Node1","Node2","Node3"]:
    k=f"{node}_rows"
    if q.get(k) != N:
        raise SystemExit(
            f"ERROR: {k}={q.get(k)} but alignment N={N}"
        )

if q["strict_contiguous_5bp_windows"] <= 0:
    raise SystemExit("ERROR: no strict 5-mer windows")

if q["coordinate_rows_ge_0.95_any_lineage"] <= 0:
    raise SystemExit("ERROR: no primary coordinate rows")

print("PASS: permanent outputs and core QC validated")
print("strict windows:",
      q["strict_contiguous_5bp_windows"])
print("coordinate rows >=0.95:",
      q["coordinate_rows_ge_0.95_any_lineage"])
PY

echo
echo "===== 8. FILE SIZES BEFORE CLEANUP ====="

ls -lh "$STATE" "$SUMMARY" "$QC" "$SITES"

echo
echo "===== 9. REMOVE BULKY ASR STATE ====="

rm -f "$STATE"

if [[ -e "$STATE" ]]; then
    echo "ERROR: state cleanup failed"
    exit 1
fi

echo "PASS: temporary ASR state removed"

echo
echo "=================================================="
echo "SUCCESS: $CHR"
echo "Finished: $(date)"
echo "=================================================="
