# gnomAD Frequencies - Usage Guide

## Overview

Query gnomAD (Genome Aggregation Database) to retrieve population allele frequencies for assessing variant rarity in rare disease analysis.

## Prerequisites

```bash
pip install requests myvariant pandas
```

## Quick Start

Tell your AI agent:
- "Get gnomAD frequency for this variant"
- "Filter my variants to keep only those with AF < 0.01"
- "Check if this variant is present in gnomAD"
- "Get population-specific frequencies from gnomAD"

## Example Prompts

### Single Variant Queries

> "What is the gnomAD allele frequency for chr7-140453136-A-T?"

> "Is rs121913527 rare in the general population?"

> "Check gnomAD exome and genome frequencies for BRAF V600E"

### Population-Specific

> "Get gnomAD frequencies by ancestry for this variant"

> "What is the East Asian frequency for this variant in gnomAD?"

> "Compare European vs African frequencies for my variants"

### Filtering

> "Filter my VCF to variants with gnomAD AF < 0.001"

> "Keep only variants absent from gnomAD"

> "Find variants rare in Europeans but common in Africans"

## What the Agent Will Do

1. Format variant identifier for gnomAD query
2. Query gnomAD API or myvariant.info
3. Extract exome and genome frequencies
4. Optionally retrieve population-specific frequencies
5. Apply filtering thresholds if requested

## Tips

- **Use both exome and genome** - take the max AF for conservative filtering
- **gnomAD v4** is the latest release with most samples (~800K exomes)
- **Population stratification** - some variants are common in one population but rare globally
- **Absent != pathogenic** - many benign variants are also rare
- **ACMG PM2** uses AF < 0.01 (1%) as supporting evidence for pathogenicity
