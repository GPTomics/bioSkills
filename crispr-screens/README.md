# crispr-screens

## Overview

Analysis of pooled CRISPR knockout and activation screens for gene essentiality and functional genomics.

**Tool type:** cli | **Primary tools:** MAGeCK, JACKS, CRISPResso2, BAGEL2

## Skills

| Skill | Description |
|-------|-------------|
| mageck-analysis | MAGeCK workflow for CRISPR screen analysis |
| jacks-analysis | JACKS for joint sgRNA efficacy and gene essentiality modeling |
| crispresso-editing | CRISPResso2 for CRISPR editing analysis |
| base-editing-analysis | Analyze base editing and prime editing outcomes |
| screen-qc | Quality control for pooled CRISPR screens |
| hit-calling | Statistical methods for identifying screen hits |
| library-design | sgRNA selection and library composition |
| batch-correction | Batch effect correction for multi-batch screens |

## Example Prompts

- "Analyze my CRISPR knockout screen with MAGeCK"
- "Compare treatment vs control screen results"
- "Assess editing efficiency with CRISPResso2"
- "Analyze my base editing experiment for C-to-T conversion"
- "Quantify prime editing outcomes from my amplicon data"
- "Identify essential genes from my dropout screen"
- "Calculate gene-level fitness scores with BAGEL2"
- "Run JACKS to model sgRNA efficacy across screens"
- "Jointly analyze multiple CRISPR screens with JACKS"

## Requirements

```bash
# MAGeCK
pip install mageck

# CRISPResso2
pip install CRISPResso2

# BAGEL2
pip install bagel

# JACKS
pip install jacks

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
