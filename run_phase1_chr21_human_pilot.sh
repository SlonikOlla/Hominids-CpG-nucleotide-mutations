#!/bin/bash
#SBATCH --job-name=hominid_chr21
#SBATCH --account=def-vasilisa
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/logs/hominid_chr21_%j.out
#SBATCH --error=results/logs/hominid_chr21_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

echo "===== START ====="
date
hostname
python3 --version

echo
echo "===== RUNNING chr21 HUMAN-LINEAGE PILOT ====="
python3 scripts/phase1_chr21_human_pilot.py

echo
echo "===== SUMMARY ====="
column -t -s $'\t' results/phase1/chr21_human_context_summary.tsv

echo
echo "===== QC ====="
column -t -s $'\t' results/phase1/chr21_human_pilot_QC.tsv

echo
echo "===== DERIVED SITE COUNT ====="
zcat results/phase1/chr21_human_derived_sites.tsv.gz | wc -l

echo
echo "===== FINISH ====="
date
