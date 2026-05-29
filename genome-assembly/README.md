# genome-assembly

## Overview

Assemble genomes from short-read, long-read, and PacBio HiFi data, then polish, scaffold, decontaminate, and quality-assess the result. Covers pre-assembly genome profiling, single-isolate through chromosome-scale and metagenome assembly, and the three-axis (contiguity/completeness/correctness) quality standard.

**Tool type:** cli | **Primary tools:** GenomeScope2, SPAdes, Flye, hifiasm, metaFlye, Pilon, YaHS, CheckM2, FCS-GX, QUAST, BUSCO, Merqury

## Skills

| Skill | Description |
|-------|-------------|
| genome-profiling | Estimate genome size, heterozygosity, ploidy, and repeat content from a k-mer spectrum (GenomeScope2/Smudgeplot) before assembling |
| short-read-assembly | De novo assembly from Illumina reads with SPAdes/MEGAHIT; the repeat-resolution ceiling and when to switch to long reads |
| long-read-assembly | Noisy ONT/PacBio CLR assembly with Flye and Canu; matching the assembler flag to the basecaller era |
| hifi-assembly | Phased, haplotype-resolved (and T2T) assembly from PacBio HiFi with hifiasm and verkko |
| metagenome-assembly | Community assembly and MAG recovery with metaFlye/metaSPAdes/MEGAHIT plus multi-binner consolidation |
| assembly-polishing | Read-type-matched consensus polishing with Racon, medaka, Polypolish, and Pilon; when not to polish |
| scaffolding | Hi-C/Omni-C and optical-map scaffolding with YaHS, SALSA2, and 3D-DNA, plus contact-map curation |
| contamination-detection | Foreign-sequence screening (FCS-GX/BlobToolKit) and MAG quality (CheckM2/GUNC/GTDB-Tk) |
| assembly-qc | Three-axis quality assessment: contiguity (QUAST/auN), completeness (BUSCO/compleasm), correctness (Merqury QV) |

## Example Prompts

- "Profile my reads with GenomeScope2 to estimate genome size and heterozygosity before assembling"
- "Assemble my bacterial genome from Illumina reads"
- "Assemble my Nanopore R10 reads with Flye"
- "Build a phased diploid assembly from HiFi reads with hifiasm"
- "Create a phased assembly with Hi-C data"
- "Assemble a metagenome and recover high-quality MAGs"
- "Polish my Nanopore assembly, but check whether the HiFi one needs polishing at all"
- "Scaffold my contigs to chromosome level with Hi-C"
- "Check my assembly for foreign contamination before submission"
- "Run QUAST, BUSCO, and Merqury QV to judge whether my assembly is good enough to publish"

## Requirements

```bash
# Profiling
conda install -c bioconda genomescope2 kmc jellyfish smudgeplot merqury

# Assemblers
conda install -c bioconda spades megahit flye canu hifiasm

# Polishing
conda install -c bioconda racon medaka polypolish pilon nextpolish

# QC
conda install -c bioconda quast busco compleasm merqury

# Scaffolding
conda install -c bioconda yahs salsa2

# Contamination / MAG quality
conda install -c bioconda checkm2 gunc gtdbtk blobtoolkit
# NCBI FCS-GX is distributed separately (see https://github.com/ncbi/fcs)
```

## Related Skills

- **genome-profiling** - First step: k-mer profiling sets every downstream expectation
- **read-qc** - Preprocess and QC reads before assembly
- **long-read-sequencing** - Basecalling, long-read QC, and alignment
- **genome-annotation** - Annotate the finished assembly
- **sequence-io** - Work with assembled FASTA files
- **workflows/genome-assembly-pipeline** - End-to-end profile -> assemble -> polish -> decontaminate -> scaffold -> QC
