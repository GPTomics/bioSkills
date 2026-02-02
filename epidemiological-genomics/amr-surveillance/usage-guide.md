# AMR Surveillance - Usage Guide

## Overview

Detect and track antimicrobial resistance genes using AMRFinderPlus for surveillance and clinical interpretation.

## Prerequisites

```bash
conda install -c bioconda ncbi-amrfinderplus
amrfinder -u  # Update database
```

## Quick Start

Tell your AI agent what you want to do:
- "Screen my genomes for AMR genes"
- "Track AMR trends across my surveillance samples"
- "Generate an AMR surveillance report"

## Example Prompts

### Gene Detection

> "Find all resistance genes in these bacterial genomes"

> "Check for carbapenemase genes in my Klebsiella isolates"

### Surveillance Analysis

> "How has ESBL prevalence changed over the past year?"

> "Are there any emerging resistance patterns in my data?"

### Clinical Interpretation

> "What antibiotics should be avoided for this isolate?"

> "Generate a resistance report for clinical review"

## What the Agent Will Do

1. Run AMRFinderPlus on genome assemblies
2. Parse results and filter high-confidence hits
3. Summarize by drug class
4. Track trends over time if longitudinal data
5. Flag emerging or critical resistance
6. Generate interpretable reports

## Tips

- **Organism flag** - Use --organism for point mutation detection
- **Quality** - Filter by coverage >90% and identity >90%
- **Database** - Update regularly with `amrfinder -u`
- **Critical genes** - Prioritize carbapenemases and MCR
- **Context** - Combine with strain typing for outbreak analysis
