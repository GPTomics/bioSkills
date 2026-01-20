# Long-Read Assembly - Usage Guide

## Overview

Long-read assemblers produce highly contiguous assemblies from ONT or PacBio reads, often achieving complete bacterial chromosomes in a single contig.

## Prerequisites

```bash
conda install -c bioconda flye canu
```

## Quick Start

Tell your AI agent what you want to do:
- "Assemble my ONT reads with Flye"
- "Run Canu assembly on my PacBio CLR data"
- "Create a complete bacterial genome from nanopore sequencing"

## Example Prompts

### ONT Assembly
> "Assemble my nanopore reads with Flye for a 5 Mb bacterial genome"
> "Run Flye on my high-coverage ONT dataset"

### PacBio Assembly
> "Assemble my PacBio CLR reads with Canu"
> "Use Flye for my PacBio data with estimated genome size of 3 Gb"

### Metagenome Assembly
> "Run Flye in metagenome mode for my ONT community sample"

## What the Agent Will Do

1. Identify input read type (ONT R9/R10, PacBio CLR)
2. Select appropriate assembler (Flye or Canu)
3. Estimate or confirm genome size
4. Run assembly with appropriate parameters
5. Convert GFA output to FASTA if needed
6. Recommend polishing steps for the assembly

## Tips

- Flye is faster and uses less memory than Canu for most applications
- Specify genome size with `--genome-size` (e.g., 5m for bacteria, 3g for mammals)
- For ONT, use `--nano-raw` for older R9 data or `--nano-hq` for R10/Q20+ data
- Long-read assemblies typically require polishing with medaka or Pilon
- Coverage of 30-50x is recommended; 100x+ for complex/repetitive genomes
- For PacBio HiFi reads, use the dedicated hifi-assembly skill instead
