# Phase 1 original Nibi workflow

## Input prerequisites

1. UCSC hg38-centered 30-primate Multiz MAF files, named `chr1.maf.gz` … `chr22.maf.gz`, staged under the original scratch `maf/` directory. Check the actual UCSC source URLs and release before running `download_autosomal_mafs.sh`.
2. The exact original seven-species fixed topology in `data/metadata/hominid_7species.tree` (not included here). **Do not invent or substitute a tree** without documenting a methodological change.
3. Nibi Slurm, `iq-tree/2.3.6`, Python 3, and suitable storage. Optional analysis/plotting requirements are listed in the root requirements files.
4. Restore the scripts to the original project's flat `scripts/` directory, or edit the paths in each driver. This code-only repository organizes files by role, while original drivers assume flat paths.

## Original run sequence

- `download_autosomal_mafs.sh`: obtain alignment inputs (inspect paths and URL).
- `make_7species_alignment.py chrN` and `make_production_map.py chrN`: create matching seven-species alignment and coordinate map.
- `run_phase1_chromosome.sh chrN`: original Slurm production driver; loads IQ-TREE 2.3.6 and calls the preceding preparation, ancestral reconstruction and `compact_asr_production.py`.
- `aggregate_phase1_genome.py`: combine 22 chromosome summary tables.
- `phase1_chromosome_validation.py`, `phase1_5mer_validation.py`, `phase1_final_audit.py`: validation. Inspect CLI usage before invoking each script.
- `phase1_motif_correlations.py`: motif-ranking conservation.
- `phase1_lineage_motif_deviations.py`: exploratory motif effect comparisons.
- `phase1_motif_multiple_testing.py`: 10,000 paired chromosome resamples; simultaneous max-T intervals over 1,536 contrasts.
- `phase1_3mer_vs_5mer.py`: 22-fold leave-one-autosome-out predictive evaluation.
- `phase1_predictive_robustness.py`: normalized log-likelihood gains and TTCTT exclusion without retraining.
- `phase1_make_figures.py`, `phase1_refine_figures.py`: plots; optional, require matplotlib and python-pptx.

## Original output locations

Production summaries: `results/phase1_genome/chrN_formal_5mer_STRICT_summary.production.tsv`; quality control and site-level output alongside these. Aggregate outputs: `results/phase1_genome/aggregated/`. Strengthening outputs: `results/phase1_strengthening/`. **None are included in this repository.**

## Validation status and caveats

The user-reported Nibi workflow passed a 22-chromosome numerical audit, but this code-only export was only statically checked during packaging; full computational reproduction was not rerun. Genome-wide data are ascertained through the hg38-centered complete-case alignment; reconstructed ancestral methylation, calibrated substitution rates, and CpG gains/losses are not part of Phase 1.
