# Short-Read Assembly - Usage Guide

## Overview

SPAdes assembles genomes from Illumina short reads using a multi-k-mer de Bruijn graph approach. It is the standard tool for bacterial, fungal, and small genome assembly.

## Prerequisites

```bash
conda install -c bioconda spades
```

## Quick Start

Tell your AI agent what you want to do:
- "Assemble my bacterial genome from paired-end Illumina reads"
- "Run SPAdes in isolate mode for a pure culture sample"
- "Extract plasmids from my bacterial sequencing data"

## Example Prompts

### Basic Assembly
> "Assemble my paired-end reads R1.fastq.gz and R2.fastq.gz with SPAdes"

### Mode-Specific Assembly
> "Run SPAdes in careful mode to minimize misassemblies"
> "Use SPAdes metagenome mode for my environmental sample"
> "Extract plasmid sequences from my bacterial Illumina data"

### Hybrid Assembly
> "Combine my Illumina and ONT reads in a hybrid SPAdes assembly"

## What the Agent Will Do

1. Verify input files exist and are paired correctly
2. Select appropriate SPAdes mode based on sample type
3. Run assembly with recommended parameters
4. Check assembly statistics (N50, number of contigs)
5. Suggest quality assessment with QUAST/BUSCO

## Tips

- Use `--isolate` mode for single bacterial isolates with uniform coverage
- Use `--careful` to reduce misassemblies at the cost of longer runtime
- SPAdes automatically selects k-mer sizes but you can specify with `-k`
- Output scaffolds.fasta is the primary result; contigs.fasta has unscaffolded sequences
- For metagenomes, prefer `--meta` mode or consider dedicated metagenome assemblers
