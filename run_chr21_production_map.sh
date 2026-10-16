#!/bin/bash
#SBATCH --job-name=chr21_prodmap
#SBATCH --account=def-vasilisa
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/logs/chr21_prodmap_%j.out
#SBATCH --error=results/logs/chr21_prodmap_%j.err

set -euo pipefail

cd /project/def-vasilisa/hominid_methyl_evolution

python3 scripts/make_chr21_production_map.py

echo
echo "===== OUTPUT ====="
ls -lh results/phase1/chr21_7species_production_map.tsv.gz

echo
echo "===== FIRST ROWS ====="
gzip -cd results/phase1/chr21_7species_production_map.tsv.gz | head -6
