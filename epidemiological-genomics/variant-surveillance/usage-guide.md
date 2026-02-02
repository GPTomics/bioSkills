# Variant Surveillance - Usage Guide

## Overview

Assign viral lineages and track variant dynamics using Nextclade for genomic surveillance of pathogens.

## Prerequisites

```bash
# Nextclade
npm install -g @nextstrain/nextclade

# Or pangolin for SARS-CoV-2
pip install pangolin
pangolin --update
```

## Quick Start

Tell your AI agent what you want to do:
- "Assign lineages to my SARS-CoV-2 sequences"
- "Track variant prevalence over time"
- "Identify emerging mutations in my surveillance data"

## Example Prompts

### Lineage Assignment

> "Run Nextclade on these viral sequences"

> "What variants are in my sequencing batch?"

### Surveillance Analysis

> "Generate a weekly variant surveillance report"

> "How has Omicron prevalence changed over the past month?"

### Mutation Tracking

> "Find the most common spike mutations in my data"

> "Are there any emerging mutations increasing in frequency?"

## What the Agent Will Do

1. Download/update Nextclade dataset
2. Run lineage assignment on sequences
3. Filter by QC status
4. Summarize lineage distribution
5. Track temporal trends if dates provided
6. Flag variants of concern and emerging mutations

## Tips

- **QC status** - Filter by 'good' or 'mediocre' QC; reject 'bad'
- **Dataset updates** - Update regularly for current lineages
- **Pangolin vs Nextclade** - Pangolin specific to SARS-CoV-2; Nextclade broader
- **Temporal analysis** - Need collection dates for trend analysis
- **VOC definitions** - WHO classifications may update
