#!/bin/bash
#SBATCH --job-name=chr22_prep
#SBATCH --account=def-vasilisa
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_prep_%j.out
#SBATCH --error=/scratch/vasilisa/hominid_methyl_evolution/logs/chr22_prep_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

echo "===== ALIGNMENT ====="
python3 scripts/make_7species_alignment.py chr22

echo
echo "===== PRODUCTION MAP ====="
python3 scripts/make_production_map.py chr22

PHY=/scratch/vasilisa/hominid_methyl_evolution/work/chr22/chr22_7species.phy
MAP=/scratch/vasilisa/hominid_methyl_evolution/work/chr22/chr22_7species_production_map.tsv.gz

echo
echo "===== CROSS-CHECK ====="

PHY_N=$(head -1 "$PHY" | awk '{print $2}')
MAP_N=$(gzip -cd "$MAP" | awk 'END{print NR-1}')

echo "PHYLIP sites: $PHY_N"
echo "Map sites:    $MAP_N"

if [[ "$PHY_N" != "$MAP_N" ]]; then
    echo "ERROR: alignment/map site-count mismatch"
    exit 1
fi

echo "PASS: alignment and production map contain identical site counts"

echo
echo "===== FILES ====="
ls -lh "$PHY" "$MAP"

echo
echo "===== DONE ====="
date
