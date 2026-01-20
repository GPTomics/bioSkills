# HiFi Assembly - Usage Guide

## Overview

hifiasm is the leading assembler for PacBio HiFi reads, producing highly contiguous and accurate genome assemblies with built-in phasing support.

## Prerequisites

```bash
# Install hifiasm
conda install -c bioconda hifiasm

# For trio binning
conda install -c bioconda yak

# For quality assessment
conda install -c bioconda quast busco seqkit
```

## Quick Start

Tell your AI agent what you want to do:
- "Assemble this HiFi dataset with hifiasm"
- "Create a phased assembly using Hi-C data"
- "Run trio-binned assembly with parental reads"

## Example Prompts

### Basic Assembly
> "Assemble this HiFi dataset with hifiasm"
> "Run hifiasm on my PacBio HiFi reads"

### Phased Assembly
> "Create a phased assembly using Hi-C data"
> "Run hifiasm with Hi-C phasing for my diploid sample"

### Trio Binning
> "Run trio-binned assembly with parental reads"
> "Use yak to create k-mer databases for trio phasing"

### Quality Assessment
> "Assess my hifiasm assembly quality with BUSCO"
> "Compare the two haplotypes from my phased assembly"

## What the Agent Will Do

1. Verify HiFi input reads and quality
2. Select assembly mode (basic, Hi-C phased, or trio)
3. Run hifiasm with appropriate parameters
4. Convert GFA output to FASTA format
5. Report assembly statistics and haplotype information
6. Suggest quality assessment and next steps

## Tips

- HiFi reads are high-accuracy; assemblies typically need no polishing
- 30-60x coverage is recommended for mammalian genomes
- Hi-C phasing produces two haplotype assemblies plus primary/alternate
- Trio binning with parental data gives the best phasing results
- Use `--n-hap` for polyploid genomes (default is 2 for diploids)
- Output GFA files can be converted to FASTA with awk or seqkit
