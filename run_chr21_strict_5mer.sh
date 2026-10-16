#!/bin/bash
#SBATCH --job-name=chr21_strict
#SBATCH --account=def-vasilisa
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/logs/chr21_strict_%j.out
#SBATCH --error=results/logs/chr21_strict_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

python3 scripts/compact_chr21_asr.py

echo
echo "===== FILES ====="
ls -lh \
 results/phase1/chr21_formal_5mer_STRICT_summary.tsv \
 results/phase1/chr21_formal_5mer_STRICT_QC.tsv

echo
echo "===== QC ====="
cat results/phase1/chr21_formal_5mer_STRICT_QC.tsv

echo
echo "===== DONE ====="
date
