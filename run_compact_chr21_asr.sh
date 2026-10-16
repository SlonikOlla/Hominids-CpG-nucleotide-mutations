#!/bin/bash
#SBATCH --job-name=hominid_5mer
#SBATCH --account=def-vasilisa
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/logs/hominid_5mer_%j.out
#SBATCH --error=results/logs/hominid_5mer_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

python3 scripts/compact_chr21_asr.py

echo
echo "===== OUTPUT FILES ====="
ls -lh \
  results/phase1/chr21_formal_5mer_summary.tsv \
  results/phase1/chr21_formal_5mer_QC.tsv

echo
echo "===== QC ====="
column -t -s $'\t' \
  results/phase1/chr21_formal_5mer_QC.tsv

echo
echo "===== FINISHED ====="
date
