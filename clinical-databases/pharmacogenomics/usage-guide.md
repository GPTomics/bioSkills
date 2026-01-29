# Pharmacogenomics - Usage Guide

## Overview

Query pharmacogenomics databases (PharmGKB, CPIC) for drug-gene interactions, interpret star alleles for metabolizer status, and identify actionable variants for drug dosing.

## Prerequisites

```bash
pip install requests pandas
```

## Quick Start

Tell your AI agent:
- "What is my CYP2D6 metabolizer status for *1/*4?"
- "Get CPIC dosing guidelines for warfarin"
- "Check if my variants affect drug metabolism"
- "Find pharmacogenomic interactions for clopidogrel"

## Example Prompts

### Star Allele Interpretation

> "What is the metabolizer status for CYP2D6 *1/*10?"

> "Interpret my CYP2C19 *2/*17 diplotype for clopidogrel"

> "Calculate activity score for CYP2D6 *4/*41"

### Drug Lookups

> "What genes affect warfarin dosing?"

> "Get CPIC guidelines for codeine and CYP2D6"

> "Which drugs are affected by DPYD variants?"

### Variant Annotation

> "Are any of my VCF variants in pharmacogenes?"

> "Check if rs4244285 affects drug metabolism"

## What the Agent Will Do

1. Query PharmGKB API for gene/drug annotations
2. Map star alleles to activity scores using CPIC tables
3. Calculate metabolizer phenotype from diplotype
4. Identify drugs with actionable PGx recommendations
5. Return dosing guidance if available

## Key Genes and Drug Classes

| Gene | Primary Drugs | Metabolism Type |
|------|---------------|-----------------|
| CYP2D6 | Codeine, tamoxifen | Activation, deactivation |
| CYP2C19 | Clopidogrel, PPIs | Activation |
| CYP2C9 | Warfarin, NSAIDs | Deactivation |
| DPYD | 5-fluorouracil | Deactivation (toxicity) |
| TPMT | Azathioprine | Deactivation (toxicity) |

## Tips

- **Activity scores** combine allele function values (0, 0.5, 1.0)
- **CPIC guidelines** provide evidence-based dosing recommendations
- **HLA alleles** require specific testing, not inferred from WGS
- **Gene duplications** can increase activity score above 2
- **PharmGKB levels** 1A/1B are highest evidence
