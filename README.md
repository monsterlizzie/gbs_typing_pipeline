
# GBS Typing Pipeline

This repository contains a fully integrated pipeline for whole-genome sequencing (WGS) analysis of *Streptococcus agalactiae* (Group B Streptococcus, GBS).  
It performs **quality control (QC)** and **genomic typing** in a single Nextflow pipeline run. Output includes a consolidated `summary.csv` covering all QC and typing results.

---

## Overview

The pipeline adapts:

- QC modules from the [Global PneumoSeq (GPS) Pipeline](https://github.com/GlobalPneumoSeq/gps-pipeline)
- Typing modules from the [GBS-Typer-sanger-nf](https://github.com/sanger-bentley-group/GBS-Typer-sanger-nf)

The result is a **modular**, **portable**, and **scalable** WGS workflow tailored for GBS genomic surveillance and characterisation.

---

## Features

The GBS Typing Pipeline includes:

### 🔬 Quality Control (QC)
| Step             | Tool(s) Used        | Description |
|------------------|---------------------|-------------|
| FASTQ validation | Internal scripts     | Paired-end and filename check |
| Read QC          | `fastp`             | Quality trimming, adapter removal |
| Assembly         | `shovill` or `unicycler` | De novo assembly |
| Assembly QC      | `quast`             | Contig count and genome size |
| Mapping QC       | `bwa`, `samtools`   | Coverage, Het-SNPs |
| Taxonomy QC      | `kraken2`, `bracken` | Confirm *S. agalactiae* identity |
| Overall QC       | Integrated rules    | PASS/FAIL based on thresholds |

### 🔬 Typing 

| Typing Task                | Tools Used                          | Description                                                                 |
|---------------------------|-------------------------------------|-----------------------------------------------------------------------------|
| **Serotyping**            | `srst2`                             | Detects capsular serotypes using GBS-SBG database |
| **Resistance Typing**     | `srst2`, `freebayes`                | IIdentifies resistance genes and GBS-specific SNP variants |
| **MLST**                  | `srst2`                             | Matches sample alleles to existing MLST scheme; detects novel alleles if present |
| **Surface Protein Typing**| `srst2`                             | Detects major surface protein genes associated with GBS virulence and immune evasion |
| **PBP Typing**     | `blastn`, `blastp`, custom scripts  | Identifies penicillin-binding protein (PBP) alleles from assemblies |

---

## Output Overview

```bash
output/
├── assemblies/                                
│   └── <sample_id>.contigs.fasta              # Assembled contigs per sample
├── sample_reports/
│   └── <sample_id>_report.csv                 # Individual sample QC reports (per sample)
├── typer/                                     # Typing-related outputs
│   ├── existing_sequence_types.txt            # Sequence types assigned by MLST
│   ├── serotype_res_incidence.txt             # Summary of serotypes and resistance presence
│   ├── drug_cat_alleles_variants.txt          # Detected resistance alleles per drug category
│   ├── gbs_res_variants.txt                   # GBS-specific resistance allele calls
│   ├── surface_protein_incidence.txt          # Detected surface proteins per sample
│   ├── surface_protein_variants.txt           # Surface protein variant-level results
│   ├── resfinder_accessions.txt               # Accessions of matched resistance genes from ResFinder
│   ├── new_mlst_alleles.log                   # Novel MLST alleles not matching existing database
│   ├── gbs_typer_report.txt                   # Combined report from all typing modules
│   ├── existing_pbp_alleles.txt               # (optional) Detected known PBP alleles from assemblies
│   ├── *_PBP_new_allele.faa                   # (optional) FASTA of novel PBP protein sequences
├── summary.csv                                # Combined summary: QC + typing per sample
```

## Typing Output Descriptions

The pipeline produces detailed output files under `output/typer/`. See full descriptions in the original [GBS-Typer documentation](https://github.com/sanger-bentley-group/GBS-Typer-sanger-nf#outputs), [Report columns details](https://docs.google.com/spreadsheets/d/1R5FFvACC3a6KCKkTiluhTj492-4cCe74HcCoklqX-X0/edit?gid=0#gid=0) including:

- `serotype_res_incidence.txt`: presence/absence of serotypes and resistance genes
- `gbs_res_variants.txt`: amino acid variants in resistance genes
- `drug_cat_alleles_variants.txt`: categorised resistance profiles
- `new_mlst_alleles.log`: log of novel MLST alleles
- `existing_sequence_types.txt`: matched STs from the MLST database
- `surface_protein_incidence.txt`, `surface_protein_variants.txt`
- `gbs_typer_report.txt`: all-in-one typing summary

---

## Key Parameters

| Parameter | Description | Default |
|----------|-------------|---------|
| `--ref_genome` | Reference genome FASTA | `data/CP129876.1.fasta` |
| `--assembler` | Assembly tool (`shovill` or `unicycler`) | `shovill` |
| `--read_len` | Read length for Bracken | `150` |
| `--pbp_contig` | Input contigs for PBP typing (optional) | `"output/assemblies/*.fasta"` |

---

## Default QC Thresholds

| Metric | Value |
|--------|-------|
| *S. agalactiae* abundance | ≥ 70.0% |
| Top non-agalactiae species | ≤ 5.0% |
| Het-SNP count | < 40 |
| Reference coverage | ≥ 70.0% |
| Assembly contigs | ≤ 500 |
| Genome size | 1.4 – 2.8 Mb |
| Depth (mean read) | ≥ 20x |

---

## Typing Thresholds

| Task | Parameter | Default |
|------|-----------|---------|
| GBS Resistance | `--gbs_res_min_coverage` | 99.9 |
|                | `--gbs_res_max_divergence` | 5 |
| MLST           | `--mlst_min_coverage` | 99.999 |
|                | `--mlst_min_read_depth` | 30 |
| Surface Typing | `--surfacetyper_min_coverage` | 99 |
|                | `--surfacetyper_max_divergence` | 8 |
|                | `--surfacetyper_min_read_depth` | 30 |
| PBP Typing     | `--pbp_frac_align_threshold` | 0.5 |
|                | `--pbp_frac_identity_threshold` | 0.5 |

---

## 🏁 How to Run

Using Docker:
```bash
nextflow run main.nf -profile standard 
```

Using Singularity:
```bash
nextflow run main.nf -profile singularity 
```

On Sanger LSF farm:
```bash
nextflow run main.nf -profile lsf 
```

Add `--run_pbptyper` and `--pbp_contig "output/assemblies/*.fasta"` if you want to run PBP typing.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 22.04
- Docker or Singularity
- Internet access (to pull containers)

---

## License

GNU General Public License v3.0

---
# Citations & Credits

This pipeline incorporates code, tools, and databases from numerous open-source sources. Below is a detailed breakdown of each component and its associated license or citation.

---

## Tools & Software

### **Nextflow**
- **Citation**: Di Tommaso et al., 2017. *Nextflow enables reproducible computational workflows*. *Nat Biotechnol* 35, 316–319.
- **License**: [Apache 2.0](https://github.com/nextflow-io/nextflow/blob/master/COPYING)

### **fastp**
- **Citation**: Chen et al., 2018. *fastp: an ultra-fast all-in-one FASTQ preprocessor*. *Bioinformatics*, 34(17): i884–i890.
- **License**: [MIT](https://github.com/OpenGene/fastp/blob/master/LICENSE)

### **Shovill**
- **Citation**: Torsten Seemann. https://github.com/tseemann/shovill
- **License**: [GPL-3.0](https://github.com/tseemann/shovill/blob/master/LICENSE)

### **Unicycler**
- **Citation**: Wick et al., 2017. *Unicycler: resolving bacterial genome assemblies from short and long sequencing reads*. *PLoS Comput Biol*.
- **License**: [GPL-3.0](https://github.com/rrwick/Unicycler/blob/main/LICENSE)

### **QUAST**
- **Citation**: Mikheenko et al., 2018. *Versatile genome assembly evaluation with QUAST-LG*, *Bioinformatics* 34(13): i142–i150.
- **License**: [GPL-2.0](https://github.com/ablab/quast/blob/master/LICENSE.txt)

### **BWA**
- **Citation**: Li H. (2013) *Aligning sequence reads with BWA-MEM*. arXiv:1303.3997v2 [q-bio.GN]
- **License**: [GPL-3.0](https://github.com/lh3/bwa/blob/master/COPYING)

### **SAMtools / BCFtools**
- **Citation**: Danecek et al., 2021. *Twelve years of SAMtools and BCFtools*. *GigaScience*, 10(2): giab008.
- **License**: [MIT](https://github.com/samtools/samtools/blob/develop/LICENSE)

### **Kraken2**
- **Citation**: Wood et al., 2019. *Improved metagenomic analysis with Kraken 2*. *Genome Biol* 20, 257.
- **License**: [MIT](https://github.com/DerrickWood/kraken2/blob/master/LICENSE)

### **Bracken**
- **Citation**: Lu et al., 2017. *Bracken: estimating species abundance in metagenomics data*. *PeerJ Comput. Sci.*
- **License**: [GPL-3.0](https://github.com/jenniferlu717/Bracken/blob/master/LICENSE)

### **SRST2**
- **Citation**: Inouye et al., 2014. *SRST2: Rapid genomic surveillance for public health and hospital microbiology labs*. *Genome Med* 6(11): 90.
- **License**: [BSD](https://github.com/katholt/srst2/blob/master/LICENSE.txt)

### **FreeBayes**
- **Citation**: Garrison and Marth. *Haplotype-based variant detection from short-read sequencing* (arXiv preprint arXiv:1207.3907)
- **License**: [MIT](https://github.com/freebayes/freebayes/blob/master/LICENSE)

### **Bedtools**
- **Citation**: Quinlan and Hall, 2010. *BEDTools: a flexible suite of utilities for comparing genomic features*. *Bioinformatics* 26(6): 841–842.
- **License**: [MIT](https://github.com/arq5x/bedtools2/blob/master/LICENSE)

---

## Databases

### **GBS Serotype Database**
- **Source**: GBS-SBG repository, commit [5e26992](https://github.com/swainechen/GBS-SBG/blob/5e26992d658b6fea3d3485ea5851141042405e23/GBS-SBG.fasta)
- **Citation**: Tiruvayipati et al., 2021. *GBS-SBG - GBS Serotyping by Genome Sequencing*. *Microb Genom*. 7(12):000688. doi: [10.1099/mgen.0.000688](https://doi.org/10.1099/mgen.0.000688)

### **GBS Resistance Gene Database**
- **Source**: `GBS_Res_Gene-DB_Final.fasta` from [Ben Metcalf’s GBS Scripts Reference repo](https://github.com/BenJamesMetcalf/GBS_Scripts_Reference/tree/master/GBS_Reference_DB), commit `feeefae`
- Contains core *S. agalactiae* AMR markers

### **PBP (β-lactam resistance) Database**
- **Files**: `GBS_bLactam_1A-DB.faa`, `GBS_bLactam_2B-DB.faa`, `GBS_bLactam_2X-DB.faa`, and `GBS_bLactam_Ref.fasta`
- **Source**: Same as above (Metcalf repo, commit `feeefae`)

### **ARG-ANNOT + ResFinder + CARD Hybrid Resistance DB**
- **File**: `ARG-ANNOT.fasta`
- **Original**: From SRST2 v0.2.0 (2016)
  - [ARGannot_r3.fasta](https://raw.githubusercontent.com/katholt/srst2/master/data/ARGannot_r3.fasta)
- **Contents**: Non-redundant resistance genes from ARG-ANNOT, ResFinder, and CARD
- **Citation**: Partridge et al., *J Antimicrob Chemother* (2018)

### **ResFinder.fasta (standalone)**
- **Contents**: 731 unique sequences from multiple ResFinder `.fsa` files
- **Source**: [ResFinder DB](https://bitbucket.org/genomicepidemiology/resfinder_db/src/d5e8eacac740da817b6d0e4015563942e1983a4b/)
- Includes: aminoglycoside, fosfomycin, macrolide, phenicol, sulphonamide, tetracyclines
- **Note**: Genes such as aph(3') variants renamed for consistency

---

## Credit Acknowledgements

- QC design adapted from: [Global PneumoSeq (GPS) Pipeline](https://github.com/GlobalPneumoSeq/gps-pipeline)
- Typing modules adapted from: [GBS-Typer-sanger-nf](https://github.com/sanger-bentley-group/GBS-Typer-sanger-nf)
- Many database files and scripts originated from: [Ben Metcalf’s GBS Scripts Repository](https://github.com/BenJamesMetcalf/GBS_Scripts_Reference)

> ⚠️ Please cite these tools/databases appropriately when using the GBS Pipeline in publications.
