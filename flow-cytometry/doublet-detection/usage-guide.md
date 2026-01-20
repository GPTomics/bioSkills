# Doublet Detection - Usage Guide

## Overview
Doublets are events where two or more cells pass through the laser simultaneously. They appear as false intermediate populations and must be removed before downstream analysis.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowDensity'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Remove doublets from my flow cytometry data"
- "Gate singlets using FSC-A vs FSC-H"
- "Identify doublets in my CyTOF data"

## Example Prompts
### Standard Doublet Removal
> "Create a singlet gate using FSC-A vs FSC-H"
> "Remove doublets from all samples in my flowSet"
> "Show the doublet rate for each sample"

### CyTOF Doublet Removal
> "Gate singlets using DNA intercalator channels"
> "Remove doublets based on Event_length"
> "Create a combined doublet filter using DNA and Event_length"

### Quality Assessment
> "Show FSC-A vs FSC-H plots before and after singlet gating"
> "Calculate the percentage of doublets removed per sample"
> "Flag samples with unusually high doublet rates"

## What the Agent Will Do
1. Identify appropriate doublet detection channels (FSC-A/H for flow, DNA/Event_length for CyTOF)
2. Create singlet gate based on pulse geometry or DNA content
3. Apply gate to remove doublets
4. Calculate doublet rates and generate QC plots
5. Return cleaned data for downstream analysis

## Tips
- FSC-A vs FSC-H is the standard method for conventional flow
- Singlets show linear A vs H relationship; doublets have higher A for given H
- CyTOF: use DNA intercalator (Ir191/Ir193) or Event_length
- Expect 1-5% doublets in PBMCs, higher in tissue digests
- High doublet rates (>15%) indicate sample preparation issues

## Detection Methods

| Method | Instrument | Principle |
|--------|------------|-----------|
| FSC-A vs FSC-H | Flow | Pulse geometry (singlets are linear) |
| FSC-A vs FSC-W | Flow | Doublets have increased width |
| DNA content | CyTOF | Doublets have ~2x DNA signal |
| Event_length | CyTOF | Doublets have longer transit time |

## Expected Doublet Rates

| Sample Type | Expected Rate |
|-------------|---------------|
| PBMCs | 1-5% |
| Cell lines | 2-10% |
| Tissue digest | 5-15% |
| Sorted cells | <1% |

## References
- flowAI: doi:10.1093/bioinformatics/btw191
- flowDensity: doi:10.1093/bioinformatics/btu677
