# Metagenome Assembly - Usage Guide

## Overview

Metagenome assembly reconstructs individual genomes from complex microbial communities. Long-read sequencing (ONT, PacBio) dramatically improves contiguity and enables recovery of complete genomes.

## Prerequisites

```bash
# Assembly
conda install -c bioconda flye spades

# Binning
conda install -c bioconda metabat2 semibin checkm2

# Taxonomy
conda install -c bioconda gtdbtk

# Utilities
conda install -c bioconda minimap2 samtools seqkit
```

## Quick Start

Tell your AI agent what you want to do:
- "Assemble this ONT metagenome with Flye"
- "Bin my metagenome assembly into MAGs"
- "Assess MAG quality with CheckM2"

## Example Prompts

### Assembly
> "Assemble this ONT metagenome with Flye"
> "Run metaSPAdes on my Illumina metagenome reads"
> "Create a hybrid assembly from short and long reads"

### Binning
> "Bin my metagenome assembly into MAGs"
> "Run MetaBAT2 binning on my assembled contigs"
> "Use SemiBin2 for deep learning-based binning"

### Quality and Taxonomy
> "Find complete circular genomes in my assembly"
> "Assess MAG quality with CheckM2"
> "Classify my MAGs with GTDB-Tk"

## What the Agent Will Do

1. Select appropriate assembler based on input data type
2. Run metagenome assembly with optimized parameters
3. Map reads back to assembly for coverage calculation
4. Bin contigs into putative genomes (MAGs)
5. Assess MAG quality with CheckM2
6. Assign taxonomy with GTDB-Tk
7. Report quality-filtered MAG statistics

## Tips

- metaFlye is recommended for long reads; metaSPAdes for Illumina only
- Long reads often recover complete circular genomes directly
- Multiple binning tools (MetaBAT2 + SemiBin2) can improve recovery
- High-quality MAGs: >90% complete, <5% contamination
- GUNC can detect chimeric MAGs missed by CheckM
- Consider co-assembly of related samples to improve binning
