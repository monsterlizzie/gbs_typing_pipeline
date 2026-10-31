# Release history

## v2.0.1 — updated 31 October 2026

- Added a QC-only profile and QC summary, with resume to full typing.
- Updated usage instructions and documented typing cut-offs and PBP reporting.
- Includes the PBS execution profile. Typing thresholds and databases are unchanged.

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

