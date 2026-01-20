# Scaffolding - Usage Guide

## Overview

Scaffolding orders and orients contigs into chromosome-level assemblies using Hi-C proximity ligation data to infer long-range contacts.

## Prerequisites

```bash
# YaHS
conda install -c bioconda yahs

# 3D-DNA
git clone https://github.com/aidenlab/3d-dna.git

# SALSA2
conda install -c bioconda salsa2

# Juicer tools
wget https://github.com/aidenlab/juicer/releases/download/v2.0/juicer_tools.2.0.jar
```

## Quick Start

Tell your AI agent what you want to do:
- "Scaffold my draft assembly using Hi-C data with YaHS"
- "Generate a contact map for Juicebox visualization"
- "Create chromosome-level scaffolds from my contigs"

## Example Prompts

### Hi-C Scaffolding
> "Scaffold my draft assembly using Hi-C data with YaHS"
> "Run 3D-DNA scaffolding with my Hi-C alignments"
> "Use SALSA2 to scaffold my contigs with proximity ligation data"

### Visualization
> "Generate a contact map for Juicebox visualization"
> "Create a .hic file from my scaffolded assembly"

### Post-Scaffolding
> "Fill gaps in scaffolds using long reads"
> "Validate my scaffolded assembly with BUSCO"
> "Check for telomeric sequences at scaffold ends"

## What the Agent Will Do

1. Align Hi-C reads to draft assembly
2. Filter contacts and remove PCR duplicates
3. Run scaffolding tool (YaHS/3D-DNA/SALSA2)
4. Generate contact map for visualization
5. Provide statistics on scaffold N50 and chromosome count
6. Suggest manual curation if needed

## Tips

- YaHS is recommended for most projects due to speed and accuracy
- Hi-C library quality is critical for scaffolding success
- Use Juicebox for manual inspection and curation of problem regions
- Telomere detection can validate chromosome completeness
- Gap-filling with long reads improves final assembly quality
- Draft assembly contiguity affects final chromosome assembly
