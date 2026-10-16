#!/bin/bash
#SBATCH --job-name=hominid_multi
#SBATCH --account=def-vasilisa
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/logs/hominid_multi_%j.out
#SBATCH --error=results/logs/hominid_multi_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

python3 scripts/phase1_chr21_multilineage.py

echo
echo "===== MULTILINEAGE RESULTS ====="
column -t -s $'\t' results/phase1/chr21_multilineage_summary.tsv

echo
echo "===== QC ====="
column -t -s $'\t' results/phase1/chr21_multilineage_QC.tsv
