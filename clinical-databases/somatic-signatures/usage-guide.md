# Somatic Mutational Signatures - Usage Guide

## Overview

Extract and analyze mutational signatures from somatic variants to identify mutagenic processes (UV, smoking, MMR deficiency, BRCA defects) using SigProfiler or MutationalPatterns.

## Prerequisites

```bash
# SigProfiler
pip install SigProfilerMatrixGenerator SigProfilerExtractor SigProfilerAssignment
```

```r
# MutationalPatterns
BiocManager::install(c('MutationalPatterns', 'BSgenome.Hsapiens.UCSC.hg38'))
```

## Quick Start

Tell your AI agent:
- "Extract mutational signatures from my tumor VCFs"
- "Fit my samples to COSMIC signatures"
- "What mutagenic processes are active in this tumor?"
- "Check for BRCA deficiency signature (SBS3)"

## Example Prompts

### Signature Extraction

> "Run SigProfiler on my VCF directory to extract de novo signatures"

> "Fit my mutation matrix to COSMIC v3 signatures"

### Interpretation

> "What signatures are dominant in this sample?"

> "Is there evidence of MMR deficiency from the signature profile?"

> "Check for UV damage signatures in my melanoma samples"

### Visualization

> "Plot the 96-context mutation spectrum for my samples"

> "Create a heatmap of signature contributions across samples"

## What the Agent Will Do

1. Generate 96-context trinucleotide mutation matrix from VCFs
2. Either extract de novo signatures (NMF) or fit to COSMIC
3. Calculate signature contributions per sample
4. Identify dominant signatures and their etiologies
5. Highlight clinical implications (HR deficiency, MMR, etc.)

## Key Signatures

| Signature | Etiology | Clinical Relevance |
|-----------|----------|-------------------|
| SBS1/SBS5 | Aging | Clock-like, correlates with age |
| SBS2/SBS13 | APOBEC | High TMB, common in breast/bladder |
| SBS3 | HR deficiency | PARP inhibitor response |
| SBS4 | Smoking | Lung cancer |
| SBS6/15/26 | MMR deficiency | Immunotherapy response |
| SBS7a/b | UV damage | Melanoma, skin cancers |

## Tips

- **WES vs WGS**: Set `exome=True` in SigProfiler for WES data
- **Minimum mutations**: Need ~50 mutations for reliable extraction
- **Cosine similarity >0.8**: Indicates same signature
- **SBS3 + high TMB**: Suggests HR deficiency with good prognosis
- **COSMIC v3.2**: Latest signature database with 79 SBS signatures
