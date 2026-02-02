# Gene Essentiality Analysis - Usage Guide

## Overview

Predict essential genes and synthetic lethal pairs through in silico gene knockout analysis with COBRApy.

## Prerequisites

```bash
pip install cobra pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Find essential genes in my metabolic model"
- "Identify synthetic lethal gene pairs"
- "Compare gene essentiality in aerobic vs anaerobic conditions"

## Example Prompts

### Essential Genes

> "Which genes are essential for growth on glucose?"

> "Find all essential genes in the E. coli model"

### Synthetic Lethality

> "Find synthetic lethal gene pairs in my model"

> "Which gene combinations are lethal when knocked out together?"

### Condition-Specific

> "Compare essential genes under aerobic and anaerobic conditions"

> "Find genes that are only essential on minimal media"

### Drug Targets

> "Identify potential drug targets from essential genes"

> "Find backup pathways through synthetic lethality analysis"

## What the Agent Will Do

1. Load the metabolic model
2. Run single gene deletion analysis
3. Classify genes as essential/non-essential
4. Optionally run double deletions for synthetic lethality
5. Report essential genes and synthetic lethal pairs
6. Compare across conditions if requested

## Tips

- **Threshold** - Use growth < 0.01 for essentiality (allows numerical tolerance)
- **Double deletions** - Very slow for large gene sets; subset first
- **Media conditions** - Essentiality depends heavily on growth conditions
- **Core vs conditional** - Distinguish always-essential from condition-specific
- **Validation** - Compare predictions with experimental data (Keio, etc.)
