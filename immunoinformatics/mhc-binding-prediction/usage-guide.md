# MHC Binding Prediction - Usage Guide

## Overview

Predict peptide-MHC binding affinity using neural network models to identify potential T-cell epitopes.

## Prerequisites

```bash
pip install mhcflurry
mhcflurry-downloads fetch
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict MHC binding for these peptides with HLA-A*02:01"
- "Scan this protein for potential epitopes"
- "Find strong binders for my patient's HLA type"

## Example Prompts

### Single Prediction

> "What is the binding affinity of SIINFEKL to HLA-A*02:01?"

> "Is this peptide a strong MHC binder?"

### Protein Scanning

> "Find all 9-mer epitopes in this spike protein"

> "Scan my antigen for epitopes binding common HLA-A alleles"

### Multiple Alleles

> "Predict binding for these peptides against all common HLA types"

> "Which of my patient's HLA alleles bind this peptide best?"

## What the Agent Will Do

1. Load MHCflurry prediction model
2. Accept peptide sequences and HLA alleles
3. Predict binding affinity (IC50) and percentile rank
4. Classify as strong/moderate/weak binder
5. Return ranked results

## Tips

- **Peptide length** - MHC-I: 8-11aa (most common 9aa); MHC-II: 13-25aa
- **Threshold** - IC50 <500nM or percentile <2% for binders
- **Patient-specific** - Use actual HLA typing for personalized predictions
- **Presentation score** - Includes processing; more biologically relevant
- **Population coverage** - 5-6 common alleles cover ~85% of population
