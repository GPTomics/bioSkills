# Genome Assembly Pipeline - Usage Guide

## Overview

This workflow orchestrates an end-to-end de novo genome assembly project and routes each step to the owning genome-assembly skill rather than restating it. It profiles the genome first (k-mer spectrum for size, heterozygosity, and ploidy), QCs reads, chooses an assembly path by data type (SPAdes for Illumina, Flye for noisy long reads, hifiasm for HiFi, metaFlye for communities), polishes only when needed, decontaminates, scaffolds when Hi-C is available, and finishes with three-axis QC that reports contiguity, completeness, and correctness together rather than N50 alone.

## Prerequisites

```bash
# Profiling + QC
conda install -c bioconda genomescope2 merqury meryl fastp nanoplot
# Assemblers
conda install -c bioconda spades flye hifiasm
# Polishing + mapping
conda install -c bioconda racon medaka minimap2 bwa samtools pilon
# Decontamination + QC (CheckM2 downloads a large DB; for MAG taxonomy add gtdbtk)
conda install -c bioconda checkm2 gunc blobtoolkit yahs quast busco
# FCS-GX is distributed separately by NCBI (container + ~470 GiB GX DB), not via bioconda: https://github.com/ncbi/fcs
```

## Quick Start

Tell your AI agent what you want to do:
- "Profile my reads then assemble the genome"
- "Assemble my bacterial genome from Nanopore reads and prove the quality"
- "Build a phased HiFi assembly and check it for false duplication"
- "Recover MAGs from my metagenome and gate them to MIMAG"

## Example Prompts

### Profiling and path choice
> "Run GenomeScope2 on my k-mer histogram first, tell me the genome size, heterozygosity, and ploidy, then recommend whether to assemble with SPAdes, Flye, or hifiasm given that I have R10 Nanopore reads."

> "I have a large heterozygous plant genome with ONT and Hi-C. Pick the assembly path and tell me where purging and scaffolding fit."

### Assembly and polishing
> "Assemble my ONT R10 reads with Flye, then polish with medaka using the model that matches my basecaller, and stop polishing at the Merqury QV plateau."

> "Build a hifiasm assembly from HiFi plus Hi-C, give me hap1 and hap2, and check whether the primary has uncollapsed haplotigs."

### Decontamination and QC
> "Screen my single-organism assembly with FCS-GX and a BlobToolKit blob plot, and tell me which foreign contigs to keep because they look host-integrated."

> "Run CheckM2 and GUNC on my MAGs and give me only the bins that pass MIMAG high-quality and a GUNC chimerism check."

> "Give me the three-axis QC for this assembly: auN/NG50 against the profiled size, compleasm completeness, and a Merqury QV from accurate reads. Do not lead with N50."

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Reads | FASTQ | Illumina (short), ONT/PacBio (long), or HiFi |
| k-mer histogram | meryl/text | from the reads, for profiling (step 0) |
| Hi-C reads | FASTQ | optional, for scaffolding (step 5) |
| Parental reads | FASTQ | optional, for trio phasing of HiFi |
| Expected size | Number | obtained from profiling, sets the NG50 denominator |

## What the Agent Will Do

1. **Profiling** - k-mer spectrum gives genome size, heterozygosity, and ploidy, setting expectations before assembling
2. **Read QC** - assesses and cleans reads; records platform and basecaller era as an assembly parameter
3. **Assembly** - routes to SPAdes, Flye, hifiasm, or metaFlye by data type
4. **Polishing** - applied only when needed; Merqury QV is the stop signal, not an iteration count
5. **Decontamination** - foreign screen for single organisms, CheckM2 plus GUNC for MAGs
6. **Scaffolding** - Hi-C scaffolding with mandatory contact-map curation
7. **Three-axis QC** - contiguity, completeness, and correctness reported together

## Assembly Strategy Selection

| Data type | Assembler | Polish? |
|-----------|-----------|---------|
| Illumina only, small/isolate | SPAdes `--isolate` | usually not |
| ONT noisy (R10) | Flye `--nano-hq` | medaka (matched model) |
| ONT noisy (R9) | Flye `--nano-raw`/`--nano-hq` per basecaller | medaka, maybe short-read |
| PacBio HiFi | hifiasm (phased) | rarely (already ~Q30+) |
| Community sample | metaFlye / metaSPAdes + binning | per-bin as needed |
| Large heterozygous eukaryote | long-read or HiFi, not short reads | per read type |

## Tips

- **Profile first.** Genome size, heterozygosity, and ploidy set the NG50 denominator, the expected haplotype count, and the purge decision; skipping this makes the assembler guess.
- **Match the flag to the chemistry.** `--nano-raw` on R10 reads silently collapses repeats while raising N50; `--nano-hq` is the R10 default.
- **A primary assembly is not a haplotype.** For allele-aware work, use phased hap1/hap2, not the primary mosaic.
- **Do not over-polish HiFi.** It is already accurate; polishing can lower QV. Measure with Merqury QV, not the reads polished with.
- **Contamination is two problems.** Cross-kingdom foreign sequence (FCS-GX, blob plot) is not the same question as intra-domain MAG chimerism (CheckM2 + GUNC).
- **Report the triad.** Scaffold N50 is glue, not sequence; a high N50 alone can hide a worse assembly. Lead with contiguity + completeness + correctness.

## Related Skills

- genome-assembly/genome-profiling - k-mer spectrum for size, heterozygosity, ploidy
- genome-assembly/short-read-assembly - SPAdes for Illumina-only genomes
- genome-assembly/long-read-assembly - Flye/Canu for noisy ONT/CLR
- genome-assembly/hifi-assembly - hifiasm phased HiFi assembly
- genome-assembly/metagenome-assembly - metaFlye/metaSPAdes + binning
- genome-assembly/assembly-polishing - Racon/medaka/Pilon, applied conditionally
- genome-assembly/contamination-detection - FCS-GX/BlobToolKit vs CheckM2/GUNC
- genome-assembly/scaffolding - YaHS Hi-C scaffolding and curation
- genome-assembly/assembly-qc - Three-axis QC with Merqury QV
- read-qc/fastp-workflow - Short-read QC before assembly
- long-read-sequencing/long-read-qc - Long-read QC and basecaller-era awareness
