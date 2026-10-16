#!/bin/bash
#SBATCH --job-name=hominid_asr
#SBATCH --account=def-vasilisa
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=results/logs/hominid_asr_%j.out
#SBATCH --error=results/logs/hominid_asr_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

module load iq-tree/3.1.2

iqtree3 \
  -s results/phase1/chr21_7species.phy \
  -st DNA \
  -m GTR+F+G4 \
  -te results/phase1/chr21_7species_model.treefile \
  -asr \
  -nt 4 \
  -pre results/phase1/chr21_7species_ASR

echo
echo "===== OUTPUT FILES ====="
ls -lh results/phase1/chr21_7species_ASR*

echo
echo "===== STATE FILE HEAD ====="
head -20 results/phase1/chr21_7species_ASR.state

echo
echo "===== FINISHED ====="
date
