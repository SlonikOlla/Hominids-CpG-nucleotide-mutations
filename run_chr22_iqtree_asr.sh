#!/bin/bash
#SBATCH --job-name=chr22_asr
#SBATCH --account=def-vasilisa
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_asr_%j.out
#SBATCH --error=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_asr_%j.err

set -euo pipefail

module load iq-tree/2.3.6

WORK=/scratch/vasilisa/hominid_methyl_evolution/work/chr22
cd "$WORK"

PHY=chr22_7species.phy
TREE=hominid_7species.tree
PREFIX=chr22_7species_ASR

echo "===== START ====="
date

echo
echo "===== INPUT ====="
head -1 "$PHY"
cat "$TREE"

echo
echo "===== IQ-TREE + ASR ====="

iqtree2 \
    -s "$PHY" \
    -m GTR+F+G4 \
    -te "$TREE" \
    -asr \
    -pre "$PREFIX" \
    -nt 1 \
    -redo

echo
echo "===== OUTPUTS ====="
ls -lh ${PREFIX}.*

echo
echo "===== ASR STATE HEADER ====="
head -3 ${PREFIX}.state

echo
echo "===== TREE ====="
cat ${PREFIX}.treefile

echo
echo "===== DONE ====="
date
