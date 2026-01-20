# Taxonomy Assignment - Usage Guide

## Overview

Taxonomic assignment classifies ASVs or OTUs to taxonomic ranks (Kingdom through Species) using reference databases like SILVA, GTDB, or UNITE.

## Prerequisites

```bash
# R packages
BiocManager::install(c('dada2', 'DECIPHER'))

# Download reference databases
# SILVA: https://zenodo.org/record/4587955
# GTDB: https://data.gtdb.ecogenomic.org/
# UNITE: https://unite.ut.ee/repository.php
```

## Quick Start

Tell your AI agent what you want to do:
- "Assign taxonomy to my ASVs using SILVA"
- "Classify my fungal ITS sequences with UNITE"

## Example Prompts

### Database Selection
> "Assign taxonomy to my ASV sequences using SILVA 138.1"

> "Use GTDB for taxonomy assignment on my environmental samples"

> "Classify my ITS sequences against the UNITE database"

### Method Selection
> "Use IDTAXA for more accurate taxonomy assignment"

> "Compare naive Bayes vs IDTAXA results"

### Confidence Filtering
> "Filter taxonomy assignments below 80% confidence at genus level"

> "What bootstrap threshold should I use for species-level calls?"

## What the Agent Will Do

1. Load ASV sequences and select appropriate reference database
2. Choose classification method (naive Bayes or IDTAXA)
3. Assign taxonomy with confidence scores
4. Filter low-confidence assignments
5. Format output as taxonomy table
6. Merge with ASV abundance table

## Tips

- SILVA is most comprehensive for 16S/18S general use
- GTDB has better taxonomy consistency for environmental samples
- UNITE is the gold standard for fungal ITS
- V4 region has limited species-level resolution
- Typical confidence thresholds: Genus 80-90%, Species 95%+

## Classification Methods

| Method | Pros | Cons |
|--------|------|------|
| Naive Bayes | Fast, handles novel sequences | May overclassify |
| Exact Matching | High precision | Misses novel taxa |
| IDTAXA | Balances speed and accuracy | Requires training |

## Reference Databases

| Database | Region | Best For |
|----------|--------|----------|
| SILVA 138.1 | 16S/18S | General bacteria/archaea |
| GTDB | 16S | Environmental, novel taxa |
| UNITE | ITS | Fungi |
| RDP | 16S | Historical compatibility |

## Confidence Thresholds

| Rank | Typical Threshold |
|------|-------------------|
| Phylum | 50-70% |
| Family | 70-80% |
| Genus | 80-90% |
| Species | 95%+ |
