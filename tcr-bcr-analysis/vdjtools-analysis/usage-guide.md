# VDJtools Analysis - Usage Guide

## Overview

Calculate immune repertoire diversity metrics, compare samples, and track clonal dynamics from TCR/BCR clonotype data.

## Prerequisites

```bash
# Download VDJtools
wget https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.zip
unzip vdjtools-1.2.1.zip

# Requires Java 8+
java -version
```

## Quick Start

Tell your AI agent:
- "Calculate Shannon diversity for my repertoire samples"
- "Find shared clonotypes between treatment and control"
- "Create a spectratype plot of CDR3 lengths"
- "Compare V gene usage across samples"

## Example Prompts

### Diversity Analysis

> "Calculate diversity metrics for all my repertoire samples"

> "What is the clonality (Gini coefficient) of each sample?"

> "Compare Shannon diversity between conditions"

### Sample Comparison

> "Find overlapping clonotypes between pre and post treatment"

> "Calculate pairwise Jaccard similarity for all samples"

> "Identify public clones shared across patients"

### Clonal Tracking

> "Track specific clonotypes across timepoints"

> "Which clones expanded after vaccination?"

> "Plot clonal dynamics over time"

## What the Agent Will Do

1. Prepare metadata file linking samples to conditions
2. Convert input to VDJtools format if needed
3. Run CalcDiversityStats for diversity metrics
4. Run overlap/comparison analyses as requested
5. Generate plots and summary tables

## Tips

- **Metadata file** is required for multi-sample analysis
- **Normalize samples** before comparing diversity
- **d50 metric** quickly shows oligoclonality (low d50 = few dominant clones)
- **F2 overlap** is frequency-weighted and more robust than Jaccard
- **Java memory** - increase with `-Xmx8g` for large datasets
