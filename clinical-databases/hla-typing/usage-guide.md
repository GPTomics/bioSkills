# HLA Typing - Usage Guide

## Overview

Call HLA alleles from NGS data (WGS, WES, or RNA-seq) for transplant matching, neoantigen prediction, and pharmacogenomic screening.

## Prerequisites

```bash
# OptiType (HLA Class I)
conda install -c bioconda optitype

# HLA-HD
# Download from https://www.genome.med.kyoto-u.ac.jp/HLA-HD/

# arcasHLA (RNA-seq)
pip install arcas-hla
arcasHLA reference --update
```

## Quick Start

Tell your AI agent:
- "Run HLA typing on my WES data"
- "Type HLA from my RNA-seq BAM"
- "Check if this patient has HLA-B*57:01 for abacavir"
- "Get 4-field HLA resolution from my sample"

## Example Prompts

### Basic Typing

> "Run OptiType on my WES BAM to get HLA-A, B, C"

> "Type HLA from my tumor RNA-seq"

### Pharmacogenomics

> "Check HLA alleles for carbamazepine contraindication"

> "Screen for HLA-B*57:01 before abacavir"

### Transplant

> "Get full HLA typing (Class I and II) for transplant matching"

## What the Agent Will Do

1. Extract reads mapping to HLA region (chr6:28-34Mb)
2. Run appropriate HLA typing tool based on data type
3. Parse results to standard nomenclature
4. Report alleles at requested resolution (2-field or 4-field)
5. Check against pharmacogenomic risk alleles if requested

## Tool Selection

| Data Type | Recommended Tool | Classes |
|-----------|------------------|---------|
| WES/WGS | OptiType, HLA-HD | I (OptiType) or I+II (HLA-HD) |
| RNA-seq | arcasHLA | I and II |
| Targeted panel | HLA-HD | I and II |

## Tips

- **OptiType is fastest** but only types Class I (HLA-A, B, C)
- **4-field resolution** (A*02:01) is clinical standard
- **Extract HLA region** from BAM before typing to speed up
- **RNA-seq typing** works well with arcasHLA from STAR-aligned BAMs
- **HLA-B*57:01 screening** is required before abacavir prescription
- **Population-specific alleles** affect pharmacogenomic risk
