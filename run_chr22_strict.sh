#!/bin/bash
#SBATCH --job-name=chr22_strict
#SBATCH --account=def-vasilisa
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_strict_%j.out
#SBATCH --error=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_strict_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

python3 scripts/compact_asr_generic.py chr22

echo
echo "===== OUTPUTS ====="
ls -lh \
 results/phase1_genome/chr22_formal_5mer_STRICT_summary.tsv \
 results/phase1_genome/chr22_formal_5mer_STRICT_QC.tsv

echo
echo "===== QC ====="
cat results/phase1_genome/chr22_formal_5mer_STRICT_QC.tsv

echo
echo "===== DONE ====="
date
