# crispr-screens

## Overview

Analysis of pooled CRISPR knockout and activation screens for gene essentiality and functional genomics.

**Tool type:** cli | **Primary tools:** MAGeCK, CRISPResso2, BAGEL2

## Skills

| Skill | Description |
|-------|-------------|
| mageck-analysis | MAGeCK workflow for CRISPR screen analysis |
| crispresso-editing | CRISPResso2 for CRISPR editing analysis |
| screen-qc | Quality control for pooled CRISPR screens |
| hit-calling | Statistical methods for identifying screen hits |
| library-design | sgRNA selection and library composition |
| batch-correction | Batch effect correction for multi-batch screens |

## Example Prompts

- "Analyze my CRISPR knockout screen with MAGeCK"
- "Compare treatment vs control screen results"
- "Assess editing efficiency with CRISPResso2"
- "Identify essential genes from my dropout screen"
- "Calculate gene-level fitness scores with BAGEL2"

## Requirements

```bash
# MAGeCK
pip install mageck

# CRISPResso2
pip install CRISPResso2

# BAGEL2
pip install bagel

# Python dependencies
pip install scipy>=1.8.0 pandas numpy matplotlib seaborn biopython scikit-learn

# Batch correction (ComBat)
pip install combat

# Additional tools
conda install -c bioconda drugz
```

## Related Skills

- **read-alignment** - Align screen reads to library
- **differential-expression** - Similar statistical concepts
- **pathway-analysis** - Enrichment of screen hits
