# Ancestral sequence-context evolution in hominids — Phase 1

Code-only archive of the 22-autosome, seven-primate ancestral sequence-context analysis. **No results, genomic data, figures, manuscripts, or derived tables are distributed here.**

## Scope

Human (hg38), chimpanzee (panTro5), bonobo (panPan2), gorilla (gorGor5), orangutan (ponAbe2), gibbon (nomLeu3), and rhesus macaque (rheMac8) were used in hg38-centered UCSC 30-primate Multiz alignments. Terminal substitution summaries are produced for human, chimpanzee, bonobo, and gorilla. The primary ancestral posterior threshold is 0.95; 0.90 and 0.99 are sensitivity thresholds. Strand-canonicalized ancestral 5-mers are counted under strict alignment-block and hg38-coordinate continuity criteria.

## Contents

- `scripts/production/`: alignment preparation, ancestral-state reconstruction, 5-mer counting, aggregation and quality-control scripts.
- `scripts/analysis/`: motif correlations, lineage-specific effects, paired chromosome-bootstrap simultaneous intervals, leave-one-autosome-out 3-mer/5-mer prediction and sensitivity analyses.
- `scripts/figures/`: original figure-generation source code (not figure files).
- `scripts/legacy_pilots/`: retained earlier chr21/chr22 development and validation scripts; **not** the final genome-wide production pipeline.
- `docs/WORKFLOW.md`: execution order, inputs, outputs and important portability limitations.

## Important reproducibility limitations

This is an archive of the **original Nibi scripts**, not a portable, turnkey workflow. Several scripts use hard-coded Nibi paths (`/project/def-vasilisa/...`, `/scratch/vasilisa/...`) and Slurm account directives. The fixed seven-species topology file (`data/metadata/hominid_7species.tree`) and UCSC alignment data are **not included**. Obtain and verify the original topology and the precise source alignment files before rerunning. Production uses IQ-TREE **2.3.6** (module `iq-tree/2.3.6`); the chr21 pilot used a different version. The original scripts were preserved rather than rewritten, so moving them into subdirectories requires adjusting script paths or copying them back into a flat `scripts/` directory for execution.

This archive does **not** claim independently verified end-to-end execution from a fresh checkout. No code changes to numerical algorithms have been made during packaging. See `docs/WORKFLOW.md` for reconstruction steps.

## Data availability

The processed Phase 1 ancestral-site datasets and numerical analysis outputs are deposited in **two complementary Zenodo records**:

- **Dataset Part 1 of 2:** [https://doi.org/10.5281/zenodo.23426694](https://doi.org/10.5281/zenodo.23426694)
- **Dataset Part 2 of 2:** [https://doi.org/10.5281/zenodo.23442145](https://doi.org/10.5281/zenodo.23442145)

**Both records are required** to reconstruct the complete ~13 GB archive. Download all archive parts and chunks from both records into the same directory, then follow `RECONSTRUCTION.txt` (included in both deposits) to reassemble the split parts, verify SHA256 checksums, reconstruct the original tar archive, and extract the data.

These records contain processed ancestral-site datasets and numerical outputs, **not** the original UCSC multiz MAF files or full production PHYLIP alignments. The deposited data do not establish ancestral DNA methylation measurements.

**Archived computational code:** [https://doi.org/10.5281/zenodo.23425960](https://doi.org/10.5281/zenodo.23425960).

## Citation and license

Publication citation, DOI and a license should be added only after the authors select and verify them. No license is granted by default.
