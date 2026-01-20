# Contamination Detection Usage Guide

## Overview

Contamination detection identifies foreign DNA in genome assemblies and assesses completeness for metagenome-assembled genomes (MAGs) and isolate assemblies.

## Tool Selection

| Tool | Purpose | Best For |
|------|---------|----------|
| CheckM2 | Completeness + contamination | All genomes |
| GUNC | Chimerism detection | MAGs |
| GTDB-Tk | Taxonomic classification | Taxonomy + novelty |
| BlobTools | Visual decontamination | Complex contamination |

## Quick Start Prompts

- "Check my MAGs for contamination using CheckM2"
- "Classify my genomes taxonomically with GTDB-Tk"
- "Detect chimeric assemblies with GUNC"
- "Remove contaminating contigs from my assembly"

## Quality Standards (MIMAG)

| Quality | Completeness | Contamination |
|---------|--------------|---------------|
| High | >90% | <5% |
| Medium | ≥50% | <10% |
| Low | <50% | <10% |

## Workflow

1. **Run CheckM2** - Completeness and contamination
2. **Run GUNC** - Detect chimeric genomes
3. **Run GTDB-Tk** - Taxonomic assignment
4. **Filter** - Apply quality thresholds
5. **Decontaminate** - Remove flagged contigs if needed

## Requirements

```bash
# CheckM2
pip install checkm2
checkm2 database --download

# GUNC
conda install -c bioconda gunc
gunc download_db .

# GTDB-Tk
conda install -c bioconda gtdbtk
gtdbtk download_gtdbtk_data
```

## Key Considerations

- **CheckM2** is faster and more accurate than original CheckM
- **GUNC** specifically detects inter-phylum chimerism
- **Combine tools** for comprehensive QC
- **MAG quality** affects downstream analyses significantly

## Related Skills

- **genome-assembly/assembly-qc** - BUSCO completeness
- **metagenomics/taxonomic-profiling** - Community analysis
