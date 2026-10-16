#!/bin/bash
#SBATCH --job-name=ape_mafs
#SBATCH --account=def-vasilisa
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/vasilisa/hominid_methyl_evolution/logs/download_mafs_%j.out
#SBATCH --error=/scratch/vasilisa/hominid_methyl_evolution/logs/download_mafs_%j.err

set -euo pipefail

DEST=/scratch/vasilisa/hominid_methyl_evolution/maf
BASE=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/maf

mkdir -p "$DEST"

for n in $(seq 1 22)
do
    chr="chr${n}"
    file="${chr}.maf.gz"

    echo "===== $chr ====="
    date

    # Reuse validated chr21 file rather than downloading again.
    if [[ "$chr" == "chr21" ]]; then
        cp -f \
          /project/def-vasilisa/hominid_methyl_evolution/data/maf/chr21.maf.gz \
          "$DEST/$file"
    else
        wget -c \
          -O "$DEST/$file" \
          "$BASE/$file"
    fi

    echo "Testing gzip..."
    gzip -t "$DEST/$file"

    echo "PASS: $file"
    ls -lh "$DEST/$file"
    echo
done

echo "===== ALL AUTOSOMES DOWNLOADED ====="
date
du -sh "$DEST"
