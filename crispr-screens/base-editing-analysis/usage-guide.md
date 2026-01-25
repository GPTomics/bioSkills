# Base Editing Analysis Usage Guide

## Overview

This guide covers analyzing base editing and prime editing outcomes from CRISPR experiments.

## Prerequisites

```bash
# CRISPResso2
pip install CRISPResso2

# Or via conda
conda install -c bioconda crispresso2
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze my base editing experiment for C-to-T conversion efficiency"
- "Quantify bystander edits in my ABE data"
- "Compare editing efficiency across multiple guides"
- "Assess prime editing outcomes from my amplicon sequencing"

## Example Prompts

### Base Editing

> "Analyze my CBE experiment. The target is a C-to-T conversion at position 6 in the editing window"

> "Calculate the editing efficiency and bystander rate for my adenine base editor samples"

### Prime Editing

> "Quantify prime editing outcomes including correct edits, partial edits, and indels"

> "Compare PE2 vs PE3 efficiency from my amplicon data"

### Quality Assessment

> "Check if my base editing has acceptable indel rates (should be <5%)"

> "Identify which samples have high bystander editing"

## What the Agent Will Do

1. Set up CRISPResso2 with base editor parameters
2. Configure expected conversions (C->T or A->G)
3. Run analysis on amplicon sequencing data
4. Extract editing efficiency and bystander metrics
5. Generate summary statistics and plots

## Tips

- Ensure amplicon covers the entire editing window
- Use high-quality paired-end reads for best results
- Include unedited control samples for comparison
- Low indel rates (<5%) indicate clean base editing
- Consider strand when interpreting conversion types
