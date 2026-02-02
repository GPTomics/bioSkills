# Neoantigen Prediction - Usage Guide

## Overview

Identify tumor neoantigens from somatic mutations for personalized cancer immunotherapy using pVACtools.

## Prerequisites

```bash
conda create -n pvactools python=3.8
conda activate pvactools
pip install pvactools
pvactools download_iedb_tools
```

## Quick Start

Tell your AI agent what you want to do:
- "Find neoantigens from my somatic VCF"
- "Predict vaccine targets from this tumor's mutations"
- "Prioritize neoantigens for my patient's HLA type"

## Example Prompts

### Basic Analysis

> "Run neoantigen prediction on my annotated VCF"

> "Find strong-binding mutant peptides in this tumor"

### Prioritization

> "Rank neoantigens by immunogenicity"

> "Find clonal neoantigens with high expression"

### Personalized

> "Use my patient's HLA type for neoantigen prediction"

> "Design a personalized cancer vaccine from these mutations"

## What the Agent Will Do

1. Verify VCF is VEP-annotated with amino acid changes
2. Run pVACseq with patient HLA alleles
3. Filter by binding affinity threshold
4. Calculate agretopicity scores
5. Prioritize by VAF and expression
6. Return ranked neoantigen candidates

## Tips

- **VEP annotation** - Required for amino acid change information
- **HLA typing** - Use patient-specific alleles (6 alleles typical)
- **Binding threshold** - IC50 <500nM standard; <50nM for strong binders
- **Clonality** - VAF >10% indicates clonal (present in most tumor cells)
- **Expression** - Unexpressed genes won't present neoantigens
