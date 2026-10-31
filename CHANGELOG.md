# Release history

## Unreleased — main

- Added `qc_only` profile / `--qc_only true` to skip all typing, including PBP.
- Added a QC-only summary with 16 QC/version columns, retaining invalid samples.
- Retained default full typing and its 40-column summary; documented cached
  continuation from QC-only to full analysis using `-resume`.
- Corrected empty report-channel fallbacks and removed duplicate typer-report
  collection. No typing thresholds, databases or genotype inference rules changed.
- Verified real Docker QC-only and resumed full execution on two real read pairs
  plus a synthetic corrupted pair. The v2.0.1 tag and release are unchanged.

## v2.0.1 — 29 October 2026

- Documented typing cut-offs, their relation to SRST2 defaults, grouped AMR
  reporting, overall-QC gating of assemblies and novel PBP sequence outputs.
- Added the missing resistance-depth and ResFinder parameters to the README.
- Included the PBS/Singularity execution profile already present on `main`
  since merge `75b24ab` but absent from the original `v2.0` tag.
- Preserved the analysis scripts, databases, typing thresholds and container
  versions from manuscript checkout `75b24ab`. No typing logic was changed.

## v2.0

Code tag: `v2.0` (commit `10b2579`; subsequently merged as `75b24ab`).

- Integrated read validation, preprocessing, assembly, mapping and taxonomy QC
  with genomic typing. Only overall-QC-passed samples enter typing.
- Integrated core PBP1A, PBP2B and PBP2X typing using QC-passed assemblies,
  with allele calls, diagnostic statuses and novel amino-acid FASTA outputs.
- Configured AMR and surface-protein gene reporting at 90% reference coverage,
  10% maximum divergence (90% identity equivalent), and 10× mean depth.
- Configured 10× serotyping and downstream MLST depth thresholds and 90%
  MLST coverage.
- Consolidated QC metrics, serotype, MLST, grouped resistance determinants,
  grouped surface proteins, PBP alleles and pipeline version in `summary.csv`.
- Distinguished `NT`, `Unresolved`, unavailable typing and corrupted inputs.
- Retained optional serotype/resistance, MLST and surface-protein modules;
  QC and PBP typing are mandatory in the integrated workflow.

