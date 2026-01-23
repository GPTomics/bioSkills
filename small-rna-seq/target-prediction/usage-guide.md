# miRNA Target Prediction - Usage Guide

## Overview

Predict target genes for miRNAs using sequence-based algorithms (miRanda, TargetScan) and validated target databases (miRTarBase).

## Prerequisites

```bash
# miRanda
conda install -c bioconda miranda

# Python packages
pip install pandas biopython gseapy
```

## Quick Start

Tell your AI agent:
- "Predict targets for my differentially expressed miRNAs"
- "Find validated targets in miRTarBase"
- "Run miRanda on my miRNA sequences"
- "Get consensus targets from multiple databases"

## Example Prompts

### Sequence-Based Prediction

> "Run miRanda to predict targets for hsa-miR-21-5p"

> "Find seed matches in 3' UTRs for my miRNA"

> "Predict targets using TargetScan scores"

### Database Lookups

> "Get experimentally validated targets from miRTarBase"

> "Query miRDB for high-confidence targets (score > 80)"

> "Find targets predicted by at least 2 databases"

### Functional Analysis

> "Run GO enrichment on predicted target genes"

> "What pathways are enriched in miR-21 targets?"

> "Find cancer-related targets of my DE miRNAs"

## What the Agent Will Do

1. Identify miRNAs for target prediction
2. Run sequence-based prediction (miRanda) or query databases
3. Filter predictions by score/confidence
4. Optionally find consensus across databases
5. Perform functional enrichment if requested

## Tips

- **Consensus targets** (2+ databases) are more reliable
- **miRTarBase** has experimentally validated targets
- **TargetScan context++ score** < -0.2 is high confidence
- **miRDB score > 80** indicates reliable prediction
- **Seed region** (nt 2-8) is critical for target recognition
- **3' UTRs** are the primary target sites
