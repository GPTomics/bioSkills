# Motif Analysis - Usage Guide

## Overview

Identify DNA sequence motifs enriched in ChIP-seq or ATAC-seq peaks to reveal transcription factor binding sites and regulatory elements.

## Prerequisites

```bash
# HOMER
conda install -c bioconda homer
perl /path/to/homer/configureHomer.pl -install hg38

# MEME Suite
conda install -c bioconda meme
```

## Quick Start

Tell your AI agent what you want to do:
- "Find enriched motifs in my ChIP-seq peaks"
- "Discover transcription factor binding sites in ATAC-seq peaks"
- "Identify co-binding factors from my peak regions"

## Example Prompts

### De Novo Discovery
> "Run de novo motif discovery on my ChIP-seq peaks"

> "Find novel motifs enriched in my ATAC-seq accessible regions"

### Known Motif Enrichment
> "Test which known TF motifs are enriched in my peaks"

> "Check if my ChIP-seq peaks contain the expected CTCF motif"

### Comparative Analysis
> "Compare motif enrichment between gained and lost peaks"

> "Find motifs specific to condition-specific peaks"

## What the Agent Will Do

1. Extract peak sequences from the genome using the peak BED file
2. Run de novo motif discovery (HOMER or MEME) on peak sequences
3. Scan peaks against known motif databases (JASPAR, HOCOMOCO)
4. Calculate enrichment statistics and rank motifs by significance
5. Generate motif logos and summary reports
6. Identify biologically relevant TFs based on enriched motifs

## Tips

- Use `-size 200` in HOMER to focus on the peak summit region
- P-value < 1e-10 indicates strong enrichment
- Target % should be much higher than background % for meaningful motifs
- If only repeat-like motifs appear, enable repeat masking
- Check that the expected TF motif ranks highly to validate ChIP quality
