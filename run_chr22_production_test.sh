#!/bin/bash
#SBATCH --job-name=chr22_prodtest
#SBATCH --account=def-vasilisa
#SBATCH --time=00:45:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_prodtest_%j.out
#SBATCH --error=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_prodtest_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

python3 scripts/compact_asr_production.py chr22

echo
echo "===== SUMMARY COMPARISON ====="

diff -u \
 results/phase1_genome/chr22_formal_5mer_STRICT_summary.tsv \
 results/phase1_genome/chr22_formal_5mer_STRICT_summary.production.tsv

echo "PASS: production summary is byte-for-byte identical"

echo
echo "===== QC CORE COMPARISON ====="

# Production QC has a few intentional extra metrics.
# Every metric from the validated QC must nevertheless have
# exactly the same value.

python3 - <<'PY'
import csv

old="results/phase1_genome/chr22_formal_5mer_STRICT_QC.tsv"
new="results/phase1_genome/chr22_formal_5mer_STRICT_QC.production.tsv"

def read(fn):
    with open(fn) as f:
        return {
            x["metric"]:int(x["count"])
            for x in csv.DictReader(f,delimiter="\t")
        }

a=read(old)
b=read(new)

bad=[]

for k,v in a.items():
    if b.get(k)!=v:
        bad.append((k,v,b.get(k)))

if bad:
    for x in bad:
        print("MISMATCH",*x)
    raise SystemExit(1)

print("PASS: every validated QC metric is identical")
print("Additional production QC metrics:",
      sorted(set(b)-set(a)))
PY

echo
echo "===== SITE TABLE ====="

ls -lh results/phase1_genome/chr22_ancestral_sites_ge0.95.tsv.gz

echo
echo "Rows including header:"
gzip -cd results/phase1_genome/chr22_ancestral_sites_ge0.95.tsv.gz | wc -l

echo
echo "Header:"
gzip -cd results/phase1_genome/chr22_ancestral_sites_ge0.95.tsv.gz | head -1 || true

echo
echo "===== DONE ====="
date
