# Differential m6A Analysis - Usage Guide

## Overview

Compare m6A methylation levels between conditions to identify epitranscriptomic changes.

## Prerequisites

```r
BiocManager::install(c('exomePeak2', 'QNB'))
```

## Quick Start

- "Compare m6A between treatment and control"
- "Find differential methylation sites"
- "Identify condition-specific m6A peaks"

## Example Prompts

### Differential Analysis

> "Run differential m6A analysis between WT and KO samples"

> "Find sites with increased m6A in treated cells"

### With Multiple Conditions

> "Compare m6A across 3 treatment groups"

> "Identify m6A changes during differentiation time course"

## What the Agent Will Do

1. Load IP and input BAMs for all conditions
2. Define experimental design
3. Call peaks and test differential methylation
4. Apply significance thresholds
5. Generate volcano plot and results table

## Tips

- **Paired design** - IP and input from same sample are paired
- **Batch effects** - Include batch in design if applicable
- **padj < 0.05** - Standard FDR threshold
- **|log2FC| > 1** - Require 2-fold change
- **Gene-level summary** - Aggregate site-level to gene-level for pathway analysis
