# Phase 2 Homo CpG / introgression / ancient methylation — release status

Status: release preparation, 2026-11-05. This note documents provenance and missing public artifacts; it is **not** a declaration that a Phase 2 release has been published.

## Purpose and scope

The follow-up study compares ancestral CpG retention in modern and archaic Homo genomes; CpG disruption between inferred Neanderthal/Denisovan-introgressed haplotypes, opposite homologues from the same person, and population controls; and reconstructed ancient methylation differences at regions annotated for archaic ancestry. Sequence-disruption and reconstructed-methylation comparisons are distinct estimands; the latter do **not** measure methylation of modern introgressed chromosomes or demonstrate epigenetic inheritance.

## Verified values transcribed from manuscript v17 (not independent reanalysis)

- Reconstructed ancestral CpG pairs: 13,743,967; selected five-genome common-callable CpGs: 3,004,026; retained across all five: 2,557,297 (85.13%).
- Corrected three-way matched site–haplotype observations: 29,717,417 Neanderthal-related and 765,630 Denisovan-related.
- Introgressed minus opposite-homologue disruption difference: +0.00452 percentage points for Neanderthal (95% two-way cluster-robust CI −0.00224 to +0.01128); +0.01176 for Denisovan (−0.02969 to +0.05320).
- Eligible conserved CpGs for ancient methylation: 2,445,803. The corrected typed-archaic exclusion-mask control pool increased from 333,529 to 706,555 eligible CpGs.
- Figure 5 is a 5-strategy × 6-contrast sensitivity matrix. Historical count-balanced Altai/Neanderthal +0.824 (Holm-adjusted p=0.019) and revised proximity-only Vindija/Neanderthal +0.773 (Holm-adjusted p=0.009) do not generalize across the stricter CpG-balanced matching strategies.

## Existing Phase 1 data/code records — NOT Phase 2 deposit confirmation

- Phase 1 ancestral data (part 1): https://doi.org/10.5281/zenodo.23426694
- Phase 1 ancestral data (part 2): https://doi.org/10.5281/zenodo.23442145
- Phase 1 computational software: https://doi.org/10.5281/zenodo.23425960

The record identities above are taken from this repository's README; live Zenodo record metadata and latest-version status were not independently verified in this update. The HGDP introgression source (Skov 2024) is https://doi.org/10.5281/zenodo.14136628; reconstructed ancient methylation source GEO GSE96833.

## Required Phase 2 outputs before a reproducible archival release

1. Original final scripts for `phase2_homo` genotype retention, 3-way matched decomposition, two-way cluster uncertainty, archaic-conditioned analyses, and methylation control/matching tests.
2. Genome-wide retention and matched-haplotype sufficient-statistics summaries, with input and output provenance/checksums.
3. A **full 30-row Figure 5 results table**, one row per strategy × ancient-donor/annotation contrast, including effect, interval, raw p, Holm-adjusted p, target/control blocks, and method/mask identifier; values should come from the original run, not from rounded heatmap labels.
4. Supplementary Table S4 joint-cluster audit, S8/S9 inputs, and revised Figure 4/5 plotting scripts plus explicit plotting data.
5. Reproduction instructions and dependency versions for Nibi, file manifest, SHA-256 checksums, scripts that verify retained sample/count totals and reported CIs.
6. License, creators, version number, release date, record relationships and final manuscript citation following verification.

Expected project workspace: `/project/def-vasilisa/hominid_methyl_evolution/phase2_homo`. The phase2 directory and individual file existence must be checked on Nibi. No inaccessible file is claimed to be present in this repository or Zenodo.

## Focus check

This release plan directly supports the principal comparison of Neanderthal-/Denisovan-introgressed versus matched non-introgressed (or no detected overlapping typed archaic segment) genomic regions for ancestral CpG retention/disruption and ancient methylation; Phase 1 context-evolution outputs remain upstream inputs rather than substitutes for Phase 2 reproducibility.
