# Metagenome Assembly Usage Guide

## Overview

Metagenome assembly reconstructs individual genomes from complex microbial communities. Long-read sequencing (ONT, PacBio) dramatically improves contiguity and enables recovery of complete genomes.

## Tool Selection

| Data Type | Recommended Tool |
|-----------|------------------|
| ONT reads | metaFlye |
| PacBio HiFi | metaFlye (--pacbio-hifi) |
| Illumina only | metaSPAdes |
| Hybrid (short + long) | metaFlye + Pilon polishing |

## Quick Start Prompts

- "Assemble this ONT metagenome with Flye"
- "Bin my metagenome assembly"
- "Find complete circular genomes"
- "Assess MAG quality with CheckM2"

## Requirements

```bash
# Assembly
conda install -c bioconda flye spades

# Binning
conda install -c bioconda metabat2 semibin checkm2

# Taxonomy
conda install -c bioconda gtdbtk

# Utilities
conda install -c bioconda minimap2 samtools seqkit
```

## Workflow Overview

```
Reads → Assembly → Mapping → Binning → QC → Taxonomy
         (Flye)    (minimap2) (MetaBAT2) (CheckM2) (GTDB-Tk)
```

## Quality Thresholds

| Category | Completeness | Contamination |
|----------|--------------|---------------|
| High-quality | >90% | <5% |
| Medium-quality | >50% | <10% |
| Low-quality | <50% | >10% |

## Related Skills

- **genome-assembly/contamination-detection** - Bin quality assessment
- **metagenomics/taxonomic-profiling** - Community composition
- **long-read-sequencing/read-qc** - Input data quality
