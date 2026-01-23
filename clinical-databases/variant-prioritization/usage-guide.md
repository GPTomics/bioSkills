# Variant Prioritization - Usage Guide

## Overview

Filter and prioritize variants from exome or genome sequencing to identify candidate disease-causing variants using population frequency, pathogenicity predictions, and clinical databases.

## Prerequisites

```bash
pip install myvariant pandas numpy
```

## Quick Start

Tell your AI agent:
- "Prioritize variants from my exome sequencing"
- "Filter to rare pathogenic variants"
- "Apply ACMG-style filtering criteria"
- "Find compound heterozygous candidates"

## Example Prompts

### Basic Filtering

> "Filter my variants to keep only those with gnomAD AF < 0.01"

> "Prioritize variants that are pathogenic in ClinVar"

> "Remove common variants and keep rare disease candidates"

### Multi-Evidence Prioritization

> "Score my variants using ClinVar, gnomAD, CADD, and REVEL"

> "Rank variants by pathogenicity evidence"

> "Apply ACMG criteria to prioritize my variant list"

### Inheritance-Based

> "Find compound heterozygous candidates for recessive disease"

> "Filter for autosomal dominant inheritance pattern"

> "Identify homozygous rare variants"

### Clinical Tiers

> "Assign clinical interpretation tiers to my variants"

> "Separate tier 1 (strong evidence) from tier 2 variants"

> "Create a clinical report with prioritized variants"

## What the Agent Will Do

1. Load variant data with annotations
2. Apply population frequency filters (gnomAD AF < 0.01)
3. Check ClinVar pathogenicity classifications
4. Score computational predictions (CADD, REVEL)
5. Optionally filter by inheritance pattern
6. Rank and tier variants by evidence strength

## Tips

- **Start with frequency filtering** - removes ~99% of variants
- **CADD phred > 20** suggests deleteriousness (top 1% of variants)
- **REVEL > 0.5** suggests pathogenicity for missense variants
- **Check ClinVar review status** - 2+ stars = higher confidence
- **Consider inheritance** - AD vs AR changes which variants matter
- **VUS requires manual review** - computational evidence is supporting only
