# Genome Scaffolding Usage Guide

## Overview

Scaffolding orders and orients contigs into chromosome-level assemblies using Hi-C proximity ligation data to infer long-range contacts.

## Tool Selection

| Tool | Strengths | Best For |
|------|-----------|----------|
| YaHS | Fast, accurate, low memory | Most projects |
| 3D-DNA | Juicebox integration | Manual curation |
| SALSA2 | Mature, well-validated | Standard workflows |

## Quick Start Prompts

- "Scaffold my draft assembly using Hi-C data with YaHS"
- "Generate a contact map for Juicebox visualization"
- "Fill gaps in scaffolds using long reads"
- "Validate my scaffolded assembly with BUSCO"

## Workflow

1. **Align Hi-C** - Map reads to draft assembly
2. **Filter contacts** - Remove PCR duplicates, low-quality alignments
3. **Scaffold** - Run YaHS/3D-DNA/SALSA2
4. **Visualize** - Generate contact map, check in Juicebox
5. **Curate** - Manual corrections if needed
6. **Gap-fill** - Close gaps with long reads
7. **Validate** - BUSCO, statistics, telomere check

## Requirements

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

## Key Considerations

- **Hi-C library quality** is critical for scaffolding success
- **Draft assembly quality** affects final chromosome assembly
- **Manual curation** often needed for problem regions
- **Telomere detection** validates chromosome completeness

## Related Skills

- **genome-assembly/long-read-assembly** - Input contigs
- **hi-c-analysis/hic-data-io** - Hi-C processing
