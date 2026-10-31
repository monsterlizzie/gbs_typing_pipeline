# GBS Typing Pipeline

The GBS Typing Pipeline is a Nextflow workflow for quality control and genomic
typing of paired-end raw reads from *Streptococcus agalactiae* (Group B
Streptococcus; GBS).

The workflow combines a GBS-focused QC pipeline, adapted from the Global
Pneumococcal Sequencing (GPS) workflow, with typing modules derived from the
Sanger GBS-Typer workflow. Version 2.0 runs QC and typing as one coordinated
analysis: samples must pass overall QC before entering the typing modules.

### Version 2.0.1 release

The [v2.0.1 release](https://github.com/monsterlizzie/gbs_typing_pipeline/releases/tag/v2.0.1)
documents the integrated v2.0 workflow and includes the PBS execution profile
added after the original `v2.0` tag. After cloning, use `git checkout v2.0.1`
for this release. The analysis scripts, databases and typing thresholds are
unchanged from the manuscript checkout `75b24ab`; use `git checkout 75b24ab`
when reproducing that exact checkout. The original `v2.0` tag is preserved.

The current `main` branch additionally supports QC-only execution and resuming
from QC-only to full typing (see below). These changes are not part of the
published `v2.0.1` tag; use `-r main` for the new mode.

The principal v2.0 updates are:

- AMR and surface-protein gene reporting uses 90% reference coverage and a
  90% sequence-identity equivalent (configured as 10% maximum divergence),
  with a separate minimum mean read depth of 10×.
- The consolidated summary groups resistance determinants by antibiotic
  category and surface-protein calls by protein family, and includes PBP calls.
- PBP1A, PBP2B and PBP2X typing is integrated into the workflow and receives
  assemblies only after the sample passes **overall QC**, not assembly QC alone.
  A separate PBP run is not required.
- Serotype reporting distinguishes named serotypes, `NT`, `Unresolved` and `NA`.

For a routine analysis, a user only needs to:

1. Put each pair of raw read files in `input/`.
2. Start Docker.
3. Run `nextflow run main.nf -profile standard` from the repository directory.
4. Open `output/summary.csv` when the run finishes.

## Workflow

```text
paired FASTQ files
        |
        v
file validation and read preprocessing
        |
        +--> read QC
        +--> assembly and assembly QC
        +--> reference mapping and mapping QC
        +--> taxonomic classification and taxonomy QC
        |
        v
overall QC
        |
        +--> failed or unusable samples: QC report only
        |
        `--> QC-passed samples
                 |
                 +--> PBP typing (full analysis only)
                 +--> serotype and resistance typing 
                 +--> MLST 
                 `--> surface-protein typing 
        |
        v
summary.csv
```

QC is the mandatory upstream gate for the v2.0 workflow. PBP typing is a core
module and always uses QC-passed assemblies produced during the same run.

## Analysis modules

### Quality control

| Stage | Main tools | Purpose |
|---|---|---|
| Input validation | `gzip`, internal scripts | Detect incomplete or malformed read files |
| Read preprocessing | fastp | Adapter removal and quality trimming |
| Read QC | Internal scripts | Assess usable sequence yield |
| Assembly | Shovill or Unicycler | Generate a draft genome assembly |
| Assembly QC | QUAST | Assess assembly length, fragmentation and depth |
| Reference mapping | BWA, SAMtools, BCFtools | Measure reference coverage and heterozygous SNPs |
| Taxonomy | Kraken 2, Bracken | Confirm that reads are predominantly GBS |
| Overall QC | Internal rules | Combine all QC results into a final status |

### Genomic typing

Typing is performed only for samples that pass overall QC.

| Module | Status in v2.0 | Main tools | Output |
|---|---|---|---|
| PBP typing | Mandatory core module | BLAST and custom scripts | PBP1A, PBP2B and PBP2X allele calls |
| Capsular serotyping | enabled by default | SRST2, GBS-SBG | `cps_type` |
| Antimicrobial-resistance typing | enabled with serotyping by default | SRST2, FreeBayes | Acquired genes and selected resistance-associated variants |
| MLST | enabled by default | SRST2 | Sequence type and allele profile |
| Surface-protein typing | enabled by default | SRST2 | Alpha-like proteins, pilus islands and other surface proteins |

## Requirements

- A Unix-like operating system
- [Nextflow](https://www.nextflow.io/) 25.04 or later
- One supported container runtime:
  - Docker for the `standard` profile
  - Singularity/Apptainer for the `singularity` and `lsf` profiles
- Internet access during initial setup to download container images and the
  configured Kraken 2 database
- Sufficient storage for assemblies, reference indices, containers and the
  Nextflow `work/` directory

The workflow has been developed for paired-end raw reads. It does not
accept assemblies as the primary workflow input.

## Installation

```bash
git clone https://github.com/monsterlizzie/gbs_typing_pipeline.git
cd gbs_typing_pipeline
```

Confirm that Nextflow and the selected container are available:

```bash
nextflow -version
docker --version
```

## Input reads

Place paired reads in `input/`, or provide another directory with `--reads`.
Both mates must share the same sample identifier.

Recognised examples include:

```text
sample_1.fastq.gz       sample_2.fastq.gz
sample_R1.fastq.gz      sample_R2.fastq.gz
sample_R1_001.fastq.gz  sample_R2_001.fastq.gz
```

The workflow also recognises `.fq`, `.fastq`, and their gzip-compressed forms.
Use one consistent naming convention per sample and do not place multiple
copies of the same read pair in the input directory.

## Running the pipeline

### Local execution with Docker

With reads in the default `input/` directory:

```bash
nextflow run main.nf -profile standard
```

With explicit input and output directories:

```bash
nextflow run main.nf \
  -profile standard \
  --reads /path/to/paired_reads \
  --output /path/to/results
```

The default output directory is `output/`.

### Local execution with Singularity/Apptainer

```bash
nextflow run main.nf \
  -profile singularity \
  --reads /path/to/paired_reads \
  --output /path/to/results
```

Container images are cached under `singularity_cache/` by default.

### Sanger FARM/LSF execution

```bash
nextflow run main.nf \
  -profile lsf \
  --reads /path/to/paired_reads \
  --output /path/to/results
```

The `lsf` profile uses the LSF executor and Singularity containers. Site-specific
queue, storage and resource settings may need to be adjusted elsewhere.

### PBS execution

The `pbs` profile uses Nextflow's `pbspro` executor with Singularity:

```bash
nextflow run main.nf \
  -profile pbs \
  --reads /path/to/paired_reads \
  --output /path/to/results
```

Adapt queue and resource settings to the cluster and validate the profile at
your site before a large run. The profile is retained from the v2.0 merge;
this documentation release does not constitute a new PBS execution test.

### Resuming an interrupted run

```bash
nextflow run main.nf \
  -profile standard \
  --reads /path/to/paired_reads \
  --output /path/to/results \
  -resume
```

Resume from the same launch directory and retain the original `.nextflow/`
history and `work/` directory. Nextflow reuses only tasks whose code, inputs,
parameters and cached work remain compatible.

## QC-only execution

Available on `main`. The `qc_only` profile runs QC without typing, including PBP,
and writes a QC-only `summary.csv`.

```bash
nextflow run monsterlizzie/gbs_typing_pipeline -r main \
  -profile standard,qc_only \
  --reads /path/to/paired_reads \
  --output qc_results \
  -name gbs_qc
```

To continue with full typing:

```bash
nextflow run monsterlizzie/gbs_typing_pipeline -r main \
  -profile standard \
  --qc_only false \
  --reads /path/to/paired_reads \
  --output full_results \
  -resume gbs_qc
```

Use the same launch directory, inputs and pipeline revision. Keep `.nextflow/`
and `work/` so completed QC tasks can be reused. The full summary is saved in
`full_results/summary.csv`. For a fresh full run, omit `-resume gbs_qc`.

## Module controls

The optional read-based typing modules are enabled by default. QC and PBP
typing are core stages and therefore do not have enable/disable parameters.

| Parameter | Default | Function |
|---|---:|---|
| `--run_sero_res` | `true` | Run capsular serotyping and resistance typing |
| `--run_surfacetyper` | `true` | Run surface-protein typing |
| `--run_mlst` | `true` | Run MLST |

Boolean parameters require an explicit value, for example:

```bash
--run_mlst false
```

The combined typer table is produced when serotype/resistance typing,
surface-protein typing and MLST are enabled together. The standard full run
enables all of these modules.

## Default QC thresholds

| Metric | Default acceptance rule |
|---|---:|
| *S. agalactiae* abundance | at least 70% |
| Top non-*S. agalactiae* species | no more than 5% |
| Reference coverage | at least 70% |
| Heterozygous SNP sites | fewer than 40 |
| Assembly contigs | no more than 500 |
| Assembly length | 1.4–2.8 Mb |
| Mean sequence depth | at least 20× |

The corresponding parameters are:

```text
--sagalactiae_percentage
--top_non_agalactiae_species_percentage
--ref_coverage
--het_snp_site
--contigs
--length_low
--length_high
--depth
```

Threshold changes should be recorded with the analysis and validated before
production use.

## Default typing thresholds

| Module | Parameter | Default |
|---|---|---:|
| Serotyping | `--serotyper_min_read_depth` | 10 |
| GBS resistance | `--gbs_res_min_coverage` | 90 |
| GBS resistance | `--gbs_res_max_divergence` | 10 |
| Additional resistance database (ResFinder) | `--other_res_min_coverage` | 90 |
| Additional resistance database (ResFinder) | `--other_res_max_divergence` | 10 |
| Resistance result processing | `--restyper_min_read_depth` | 10 |
| MLST | `--mlst_min_coverage` | 90 |
| MLST | `--mlst_min_read_depth` | 10 |
| Surface proteins | `--surfacetyper_min_coverage` | 90 |
| Surface proteins | `--surfacetyper_max_divergence` | 10 |
| Surface proteins | `--surfacetyper_min_read_depth` | 10 |
| PBP target detection | `--pbp_frac_align_threshold` | 0.5 |
| PBP target detection | `--pbp_frac_identity_threshold` | 0.5 |

Coverage and divergence values above are percentages; read-depth values are
mean depth in ×. The 90% coverage and 10% maximum-divergence settings match
[SRST2's default gene-reporting parameters](https://github.com/katholt/srst2/blob/master/scripts/srst2.py).
For reporting, the latter is described as a 90% sequence-identity equivalent;
the actual command-line parameter remains `--max_divergence`.
The pipeline additionally filters AMR and surface-protein calls at 10× mean
depth. This is stricter than SRST2's default 5× low-depth warning threshold:
SRST2 flags dubious calls, whereas these pipeline result processors filter
gene-presence calls. These typing thresholds are separate from the 20×
genome-wide QC depth requirement.

`--mlst_min_read_depth` controls downstream MLST allele processing; it does not
replace SRST2's own depth-warning setting. PBP typing is assembly-based and
uses its separate fractional thresholds, not the SRST2 gene-reporting cut-offs.

## Output files

The standard full run produces:

```text
output/
├── summary.csv
├── assemblies/
│   └── <sample_id>.contigs.fasta
├── qc_reports/
│   └── <sample_id>_qc.csv
├── typer/
│   ├── gbs_typer_report.txt
│   ├── serotype_res_incidence.txt
│   ├── drug_cat_alleles_variants.txt
│   ├── gbs_res_variants.txt
│   ├── resfinder_accessions.txt
│   ├── existing_sequence_types.txt
│   ├── new_mlst_alleles.log
│   ├── surface_protein_incidence.txt
│   ├── surface_protein_variants.txt
│   └── existing_pbp_alleles.txt
├── pbp_target_status.txt
├── pbp_allele_status.txt
├── <sample_id>_<PBP_target>_PBP_new_allele.faa  # when applicable
└── <sample_id>_new_mlst_alleles.fasta   # when applicable
```

`summary.csv` is the principal analysis table. It combines QC status, QC
metrics, serotype, MLST, grouped resistance determinants, grouped
surface-protein calls, PBP alleles and pipeline version.


### Summary report and detailed typing report

The pipeline provides two reports for different purposes:

| Report | Intended use | Contents |
|---|---|---|
| `output/summary.csv` | Routine review, surveillance and downstream analysis | QC results, grouped antimicrobial-resistance determinants, serotype, ST, grouped surface-protein calls, PBP1A/PBP2B/PBP2X allele calls and pipeline version |
| `output/typer/gbs_typer_report.txt` | Inspection of individual typing calls | Serotype, ST, the seven MLST loci, individual resistance-gene incidence, selected GBS-specific resistance variants, surface-protein incidence and pipeline version |

The combined detailed typing report is generated when serotype/resistance typing, MLST
and surface-protein typing are all enabled. Its individual gene columns retain
more detail than the grouped positive/negative determinant columns in
`summary.csv`. PBP results are not part of `gbs_typer_report.txt`; they are
reported in `summary.csv`, `typer/existing_pbp_alleles.txt`, the PBP status
files and any novel PBP FASTA output.

The grouped AMR columns contain detected determinant names or variants, or
`neg` when no included determinant is reported. They are **not phenotypic
susceptibility categories**: `neg` does not establish susceptibility, and a
determinant call should be interpreted using the relevant inference rules.
Unavailable results must not be treated as negative calls. The summary covers
erythromycin, clindamycin, tetracycline, chloramphenicol, high-level gentamicin
resistance, fluoroquinolones, vancomycin and other aminoglycosides.

Column descriptions for the detailed typing report are available in the
[GBS-Typer report-column reference](https://docs.google.com/spreadsheets/d/1R5FFvACC3a6KCKkTiluhTj492-4cCe74HcCoklqX-X0/edit?usp=sharing).
In that reference, use rows where the `category` column is
`in_silico_analysis`.

## Interpreting QC and serotype values

Typing is attempted only for samples with `Overall_QC=PASS`.

| Value | Interpretation |
|---|---|
| A named serotype, for example `Ia`, `Ib`, `II` or `III` | A reportable in-silico serotype was assigned |
| `NT` | SRST2 detected depth was below the configured serotyping depth threshold ( default : 10x) |
| `Unresolved` | Overall QC passed, but SRST2 did not produce a valid reportable serotype/fullgenes call |
| `NA` | Serotyping was unavailable because overall QC failed or the input reads were unusable |

Corrupted inputs are reported explicitly in `Overall_QC`, for example
`READ_ONE_CORRUPTED` or `READ_TWO_CORRUPTED`, and their unavailable typing
fields remain `NA`.

### Interpreting PBP allele values

PBP typing examines PBP1A, PBP2B and PBP2X in the QC-passed de novo assembly.
The three `pbp` columns in `summary.csv` contain either an allele number or a
diagnostic value.

| Value | Interpretation |
|---|---|
| An allele number | The translated PBP sequence exactly matched an allele in the current database |
| `NEW` | A qualifying translated PBP target fragment was detected but did not exactly match a known database allele; the sequence is also written to `<sample_id>_<PBP_target>_PBP_new_allele.faa` |
| `PARTIAL_PBP_GENE` | The predicted PBP region crosses the start or end of an assembly contig, so a complete allele cannot be assigned |
| `PBP_NOT_DETECTED` | No target passed the initial PBP detection thresholds |
| `NO_BLAST_HIT` | A complete query was produced, but the protein search returned no database hit |
| `NO_ALLELE_MATCH` | The module could not make a qualifying known- or new-allele call |
| `MODULE FAILURE` | A computational step in PBP typing failed; inspect the Nextflow log and PBP status files |
| `NA` | PBP typing was unavailable because the sample did not pass overall QC or its reads were unusable |

`PARTIAL_PBP_GENE` usually indicates that the draft assembly is fragmented at
that locus. It should not automatically be interpreted as a genuinely
disrupted gene in the isolate. Likewise, `PBP_NOT_DETECTED`, `NO_BLAST_HIT` and
`NO_ALLELE_MATCH` are diagnostic outcomes rather than allele names. Review
`pbp_target_status.txt`, `pbp_allele_status.txt` and the assembly before making
a biological interpretation.

New PBP files contain amino-acid sequences (`.faa`). The target names are
`GBS1A-1`, `GBS2B-1` and `GBS2X-1`. Preserve these files with the summary and
assemblies: a `NEW` label and contig coordinates alone do not contain the
sequence. Filenames use the sample identifier passed to the module; historical
standalone runs may include `.contigs` in that identifier. Different samples
labelled `NEW` may carry the same sequence, so the number of `NEW` calls is not
the number of distinct novel alleles.

## Storage management

Nextflow's `work/` directory can be substantially larger than the published
results. Do not remove it while a run is active or while `-resume` may still be
needed.

After a run has been validated and all required outputs have been safely
copied, cached task data can be reviewed with Nextflow's cleanup facilities.
Always confirm the target run before deleting cached work.

The `--lite true` option can be used where supported by the workflow to reduce
retained intermediate output. Validate storage-saving settings on a small run
before production use.

## Troubleshooting

### No input pairs are detected

Check that:

- both pairs files are in the directory supplied with `--reads`;
- filenames use a supported `_1/_2` or `_R1/_R2` pattern;
- sample identifiers match exactly between mates;
- duplicate compressed and uncompressed copies are not present.

### Docker processes exit immediately

Confirm that Docker Desktop/the Docker daemon is running:

```bash
docker info
```

Then rerun with `-resume` after correcting the problem.

### A read is reported as corrupted

Validate gzip-compressed reads independently:

```bash
gzip -t sample_1.fastq.gz
gzip -t sample_2.fastq.gz
```

Redownload or restore invalid files. Do not decompress and recompress a damaged
file as a substitute for recovering the original sequence data.

### A process fails

Inspect `.nextflow.log` and the process work directory reported by Nextflow.
Useful files include:

```text
.command.sh
.command.out
.command.err
.exitcode
```

## Project history and attribution

The GBS Typing Pipeline integrates two earlier bodies of work:

- QC design and modules adapted from the
  [Global Pneumococcal Sequencing pipeline](https://github.com/GlobalPneumoSeq/gps-pipeline)
- GBS typing modules adapted from
  [GBS-Typer-sanger-nf](https://github.com/sanger-bentley-group/GBS-Typer-sanger-nf)

Version 2.0 represents the integrated QC-plus-typing workflow. QC is mandatory; full analysis includes PBP typing, while the read-based typing modules are
configurable.

## Citations

Users should cite this pipeline and the underlying methods relevant to their
analysis. Principal software and database references include:

- Di Tommaso P, et al. Nextflow enables reproducible computational workflows.
  *Nature Biotechnology*. 2017;35:316–319.
- Chen S, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor.
  *Bioinformatics*. 2018;34:i884–i890.
- Wick RR, et al. Unicycler: resolving bacterial genome assemblies from short
  and long sequencing reads. *PLoS Computational Biology*. 2017.
- Mikheenko A, et al. Versatile genome assembly evaluation with QUAST-LG.
  *Bioinformatics*. 2018;34:i142–i150.
- Li H. Aligning sequence reads, clone sequences and assembly contigs with
  BWA-MEM. arXiv:1303.3997. 2013.
- Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience*.
  2021;10:giab008.
- Wood DE, et al. Improved metagenomic analysis with Kraken 2. *Genome
  Biology*. 2019;20:257.
- Lu J, et al. Bracken: estimating species abundance in metagenomics data.
  *PeerJ Computer Science*. 2017.
- Inouye M, et al. SRST2: rapid genomic surveillance for public health and
  hospital microbiology labs. *Genome Medicine*. 2014;6:90.
- Tiruvayipati S, et al. GBS-SBG—GBS serotyping by genome sequencing.
  *Microbial Genomics*. 2021;7:000688.
- Garrison E, Marth G. Haplotype-based variant detection from short-read
  sequencing. arXiv:1207.3907.
- Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing
  genomic features. *Bioinformatics*. 2010;26:841–842.

Reference data were derived from GBS-SBG, the Ben Metcalf GBS Scripts Reference
repository, ResFinder, ARG-ANNOT and associated resources. Consult the source
repositories and database files for version-specific provenance and licensing.

## License

This project is distributed under the GNU General Public License v3.0. See
[`LICENSE`](LICENSE) for the full terms.
