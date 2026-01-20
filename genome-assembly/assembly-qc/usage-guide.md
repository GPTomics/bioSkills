# Assembly QC - Usage Guide

## Overview

Assembly QC evaluates the quality of genome assemblies using contiguity metrics (QUAST) and gene completeness (BUSCO).

## Prerequisites

```bash
conda install -c bioconda quast busco
```

## Quick Start

Tell your AI agent what you want to do:
- "Check the quality of my genome assembly"
- "Run BUSCO on my bacterial assembly"
- "Compare two assemblies with QUAST"

## Example Prompts

### Basic Quality Check
> "Run QUAST on my assembly to get N50 and other statistics"
> "Check my assembly quality with BUSCO using the bacteria lineage"

### Comparative Analysis
> "Compare my SPAdes and Flye assemblies with QUAST"
> "Run QUAST with a reference genome to detect misassemblies"

### Lineage-Specific BUSCO
> "Run BUSCO on my fungal assembly using fungi_odb10"
> "Check gene completeness of my plant genome with viridiplantae lineage"

## What the Agent Will Do

1. Run QUAST to calculate contiguity statistics (N50, L50, total length)
2. Run BUSCO with appropriate lineage database
3. Report key metrics and quality assessment
4. Flag potential issues (low N50, missing BUSCOs, contamination)
5. Suggest next steps based on quality results

## Tips

- Target >95% BUSCO completeness for high-quality assemblies
- N50 should be interpreted relative to expected chromosome/genome size
- Use `busco --list-datasets` to find the appropriate lineage database
- QUAST with a reference (`-r`) can detect misassemblies and structural errors
- Run both tools together for a comprehensive quality picture
- High duplication in BUSCO may indicate contamination or assembly issues
