# scirpy Analysis - Usage Guide

## Overview

Integrate single-cell TCR/BCR data with gene expression analysis using scirpy, enabling joint analysis of immune receptor sequences and cellular phenotypes.

## Prerequisites

```bash
pip install scirpy scanpy mudata
```

## Quick Start

Tell your AI agent:
- "Load 10x VDJ data into my AnnData object"
- "Find clonally expanded cells"
- "Compare repertoire diversity between conditions"
- "Plot V gene usage by cell type"

## Example Prompts

### Data Loading

> "Add 10x VDJ annotations to my scRNA-seq data"

> "Load TCR data from filtered_contig_annotations.csv"

> "Convert my AnnData to AIRR format"

### Clonal Analysis

> "Define clonotypes by CDR3 amino acid sequence"

> "Which clonotypes are expanded (>10 cells)?"

> "Track clonotypes across timepoints"

### Integration

> "Find genes differentially expressed in expanded clones"

> "Color my UMAP by clonal expansion"

> "Compare cell types with and without TCR"

### Visualization

> "Plot V-J gene usage pairing"

> "Create a spectratype plot by cell type"

> "Show repertoire overlap between samples"

## What the Agent Will Do

1. Load VDJ data and integrate with AnnData
2. Run chain QC to identify doublets/orphans
3. Define clonotypes by CDR3 identity
4. Calculate clonal expansion
5. Compute diversity metrics
6. Generate visualizations

## Tips

- **Chain pairing QC** - filter multichain cells (likely doublets)
- **Dual IR handling** - decide how to treat cells with 2 alpha/beta chains
- **Integration** - scirpy uses same AnnData as Scanpy
- **Memory** - large datasets may need subsetting
- **AIRR export** - for compatibility with other tools
