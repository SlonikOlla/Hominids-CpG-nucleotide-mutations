#!/bin/bash
#SBATCH --job-name=hominid_tree
#SBATCH --account=def-vasilisa
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=results/logs/hominid_tree_%j.out
#SBATCH --error=results/logs/hominid_tree_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

module load iq-tree/3.1.2

iqtree3 \
  -s results/phase1/chr21_7species.phy \
  -st DNA \
  -m GTR+F+G4 \
  -te data/metadata/hominid_7species.tree \
  -nt 4 \
  -pre results/phase1/chr21_7species_model

echo
echo "===== TREE ====="
cat results/phase1/chr21_7species_model.treefile

echo
echo "===== MODEL SUMMARY ====="
grep -E \
"Model of substitution|Base frequencies|Rate parameter|Gamma shape|Tree length" \
results/phase1/chr21_7species_model.iqtree || true

echo
echo "===== FINISHED ====="
date
