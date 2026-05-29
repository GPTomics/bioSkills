# Metagenome Assembly - Usage Guide

## Overview

Metagenome assembly reconstructs individual microbial genomes (MAGs) from mixed-community sequencing. Unlike isolate assembly, the assembler cannot assume uniform coverage: abundant and rare species coexist across a wide depth range, and co-occurring strains collapse into consensus contigs. The deliverable is therefore a set of consensus population genomes (MAGs), gated against the MIMAG quality standard, not a single contiguous assembly. Modern practice assembles with a metagenome-aware assembler (metaFlye for ONT, metaSPAdes/MEGAHIT for Illumina, hifiasm-meta/metaMDBG for HiFi), runs multiple binners, consolidates with DAS_Tool, then QCs every MAG with CheckM2 and GUNC and classifies with GTDB-Tk. HiFi/long reads are the biggest quality jump available, recovering the complete circular genomes and rRNA operons short-read MAGs lose.

## Prerequisites

```bash
# Assembly
conda install -c bioconda flye spades megahit
# HiFi metagenome assemblers
pip install metaMDBG   # or conda; hifiasm-meta from bioconda

# Binning + consolidation
conda install -c bioconda metabat2 maxbin2 concoct semibin vamb das_tool

# MAG QC + taxonomy (each downloads a large reference DB)
conda install -c bioconda checkm2 gunc gtdbtk   # GTDB-Tk release must match the DB

# Strain resolution + mapping utilities
conda install -c bioconda instrain minimap2 samtools seqkit
```

## Quick Start

Tell your AI agent what you want to do:
- "Assemble this ONT metagenome with metaFlye"
- "Recover MAGs from my Illumina gut microbiome co-assembly"
- "Run multiple binners and consolidate with DAS_Tool"
- "QC my MAGs against MIMAG with CheckM2 and GUNC"
- "Resolve strains in my recovered MAGs with inStrain"

## Example Prompts

### Assembly
> "Assemble this ONT gut-microbiome dataset with metaFlye in --meta mode, then tell me whether I need to polish before binning."
> "My soil Illumina co-assembly won't fit in RAM with metaSPAdes -- switch to MEGAHIT and explain the contiguity tradeoff."
> "I have PacBio HiFi from an ocean sample; pick hifiasm-meta or metaMDBG and explain why HiFi recovers circular MAGs short reads can't."

### Binning and consolidation
> "Map all four of my samples back to the assembly, bin with MetaBAT2, SemiBin2, and MaxBin2, then consolidate with DAS_Tool."
> "I only have one sample -- explain why my MetaBAT2 bins are poor and what would actually fix it."
> "Set up differential-coverage binning across my time-series samples."

### QC, taxonomy, and strains
> "QC my DAS_Tool bins with CheckM2 and GUNC and tell me which pass MIMAG high-quality."
> "My short-read MAG is 95% complete and 2% contaminated but isn't HQ -- why, and does HiFi fix it?"
> "Classify my MAGs with GTDB-Tk and dereplicate redundant genomes across samples."
> "Use inStrain to check whether a strain-level claim from one consensus MAG holds up."

## What the Agent Will Do

1. Select a metagenome-aware assembler from read type (metaFlye for ONT, metaSPAdes/MEGAHIT for Illumina, hifiasm-meta/metaMDBG for HiFi) and run it in `--meta` mode.
2. Polish ONT/CLR contigs if needed (HiFi usually needs none).
3. Map reads back per sample to build differential-coverage depth profiles for binning.
4. Run two or more binners (composition + coverage), then consolidate with DAS_Tool into a non-redundant set.
5. QC every MAG with CheckM2 (completeness/contamination) and GUNC (chimerism).
6. Classify passing MAGs with GTDB-Tk and dereplicate across samples.
7. Resolve strains with inStrain when the question is strain-level.
8. Report MAG count by MIMAG tier and community fraction binned -- not a single N50.

## Tips

- Always use `--meta` mode; isolate assemblers treat the community's coverage spread and strain bubbles as errors and produce garbage.
- metaSPAdes accepts exactly one paired library -- concatenate libraries first or use MEGAHIT for many.
- Differential-coverage binning needs multiple samples with abundance variation; with one sample, binning collapses to weak composition-only signal. Add samples, do not retune.
- Never run a single binner: each recovers a partially-different genome set. Run several and consolidate with DAS_Tool.
- Pair CheckM2 with GUNC. A high-completeness, low-contamination MAG can still be a chimera that only GUNC catches.
- Short-read MAGs routinely meet 90/5 yet fail MIMAG HQ because the rRNA operon collapsed and stayed unbinned. HiFi/long reads span the operon and fix it.
- The rare biosphere does not assemble (below ~3-5x coverage). A MAG catalogue is a biased sample of the abundant fraction; complement with read-based profiling.
- Match the GTDB-Tk reference-package release to the installed binary, or classification silently fails.
- A MAG is a population consensus. Treat any strain-level claim from a single MAG as suspect until inStrain or strain-aware assembly confirms it.

## Related Skills

- contamination-detection - CheckM2/GUNC interpretation and chimerism forensics for MAGs
- assembly-qc - Isolate-assembly QC; the uniform-coverage assumption metagenomes abandon
- assembly-polishing - Polish ONT/CLR meta-contigs before binning
- metagenomics/kraken-classification - Read-based taxonomy; recovers the rare biosphere assembly cannot
- metagenomics/abundance-estimation - Community abundance downstream of recovered MAGs
- long-read-sequencing/long-read-qc - Read-level QC and host removal before assembly
