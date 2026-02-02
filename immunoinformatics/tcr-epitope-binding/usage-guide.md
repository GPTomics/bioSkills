# TCR-Epitope Binding - Usage Guide

## Overview

Predict TCR-epitope specificity using deep learning models and database matching to identify antigen-reactive T-cells.

## Prerequisites

```bash
pip install pandas torch scikit-learn
# ERGO-II: git clone https://github.com/IdoSpringer/ERGO-II
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict what antigens this TCR sequence recognizes"
- "Match my TCRs to known epitopes in VDJdb"
- "Cluster TCRs that likely share specificity"

## Example Prompts

### Specificity Prediction

> "What epitopes might this CDR3 beta sequence recognize?"

> "Predict binding between my TCRs and these candidate epitopes"

### Database Matching

> "Find matches for my TCRs in VDJdb"

> "Identify TCRs recognizing viral epitopes"

### Repertoire Analysis

> "What fraction of my repertoire recognizes known antigens?"

> "Cluster TCRs by predicted specificity"

## What the Agent Will Do

1. Parse TCR sequence data (CDR3 alpha/beta)
2. Match against VDJdb or other databases
3. Run ERGO-II prediction if available
4. Cluster similar TCRs
5. Report potential epitope specificities

## Tips

- **CDR3 beta** - Most informative for specificity; alpha adds ~20%
- **VDJdb** - Curated database of known TCR-epitope pairs
- **Clustering** - TCRs within 1-3 edit distance often share specificity
- **ERGO-II** - Deep learning provides better predictions than simple matching
- **Validation** - Predicted specificities should be validated experimentally
