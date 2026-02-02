# Guide RNA Design - Usage Guide

## Overview

Design and score CRISPR guide RNAs for Cas9 and Cas12a nucleases using on-target activity prediction models.

## Prerequisites

```bash
pip install biopython crisprscan
```

## Quick Start

Tell your AI agent what you want to do:
- "Design guides to knock out BRCA1"
- "Find sgRNAs targeting exon 3 of TP53 with high activity scores"
- "Design Cas12a guides for this target sequence"

## Example Prompts

### Gene Knockout

> "Design the best 3 guides to knock out the KRAS gene"

> "Find sgRNAs in the first coding exon of MYC"

### Custom Targets

> "Score this guide sequence for on-target activity: ATCGATCGATCGATCGATCG"

> "Design guides within 100bp of position chr7:140453136"

### Different Nucleases

> "Design Cas12a guides for my target region"

> "Find SaCas9-compatible guides for AAV delivery"

## What the Agent Will Do

1. Retrieve or accept the target gene/sequence
2. Identify all valid PAM sites in the region
3. Extract guide sequences for each PAM
4. Score guides using activity prediction models
5. Filter by GC content and avoid poly-T stretches
6. Return top-ranked guides with scores and positions

## Tips

- **Exon targeting** - Target early coding exons for reliable knockout
- **GC content** - Optimal range is 40-70%; outside this reduces activity
- **Avoid poly-T** - Four or more consecutive T's terminate transcription
- **Multiple guides** - Design 3-5 guides per gene for redundancy
- **Off-target check** - Always check off-targets before ordering guides
